from sympy import Symbol, solve, simplify, Matrix
import sympy
import copy
import pprint
from itertools import chain

class Circ():
    """
    All nodes and values are not case Sensitive (defaults to Capitol followed by lower case)
    Gnd, jw, kbt, Isource: reserved names

    Voltages start with V
    Inductors start with L
    Resistors start with R
    Capacitors start with C

    """
    jw = Symbol("jw")
    w = Symbol("w")
    j = Symbol("j")
    kbt = Symbol("kbt")

    def __init__(self):
        self.nodes = dict()
        self.var_list = dict()

    def add_res(self, res, plus_term, minus_term):
        res = Circ.format_name(res)
        plus_term = Circ.format_name(plus_term)
        minus_term = Circ.format_name(minus_term)

        if res not in self.var_list:
            self.var_list[res] = Symbol(res)
        if plus_term not in self.var_list:
            self.var_list[plus_term] = Symbol(plus_term)
        if minus_term not in self.var_list:
            self.var_list[minus_term] = Symbol(minus_term)

        if plus_term not in self.nodes:
            self.nodes[plus_term] = []
        if minus_term not in self.nodes:
            self.nodes[minus_term] = []

        self.nodes[plus_term] += [StoredOp(plus_term, minus_term, res, top=False, add_jw=False, neg=True)]
        self.nodes[minus_term] += [StoredOp(plus_term, minus_term, res, top=False, add_jw=False, neg=False)]

    def add_cond(self, cond, plus_term, minus_term):
        self.add_vccs(cond, plus_term, minus_term, minus_term, plus_term)

    def add_cap(self, cap, plus_term, minus_term):
        cap = Circ.format_name(cap)
        plus_term = Circ.format_name(plus_term)
        minus_term = Circ.format_name(minus_term)

        if cap not in self.var_list:
            self.var_list[cap] = Symbol(cap)
        if plus_term not in self.var_list:
            self.var_list[plus_term] = Symbol(plus_term)
        if minus_term not in self.var_list:
            self.var_list[minus_term] = Symbol(minus_term)

        if plus_term not in self.nodes:
            self.nodes[plus_term] = []
        if minus_term not in self.nodes:
            self.nodes[minus_term] = []

        self.nodes[plus_term] += [StoredOp(plus_term, minus_term, cap, top=True, add_jw=True, neg=True)]
        self.nodes[minus_term] += [StoredOp(plus_term, minus_term, cap, top=True, add_jw=True, neg=False)]

    def add_ind(self, ind, plus_term, minus_term):
        ind = Circ.format_name(ind)
        plus_term = Circ.format_name(plus_term)
        minus_term = Circ.format_name(minus_term)

        if ind not in self.var_list:
            self.var_list[ind] = Symbol(ind)
        if plus_term not in self.var_list:
            self.var_list[plus_term] = Symbol(plus_term)
        if minus_term not in self.var_list:
            self.var_list[minus_term] = Symbol(minus_term)

        if plus_term not in self.nodes:
            self.nodes[plus_term] = []
        if minus_term not in self.nodes:
            self.nodes[minus_term] = []

        self.nodes[plus_term] += [StoredOp(plus_term, minus_term, ind, top=False, add_jw=True, neg=True)]
        self.nodes[minus_term] += [StoredOp(plus_term, minus_term, ind, top=False, add_jw=True, neg=False)]

    def add_vccs(self, gain, plus_term, minus_term, ref_plus, ref_minus):
        gain = Circ.format_name(gain)
        plus_term = Circ.format_name(plus_term)
        minus_term = Circ.format_name(minus_term)
        ref_plus = Circ.format_name(ref_plus)
        ref_minus = Circ.format_name(ref_minus)

        if gain not in self.var_list:
            self.var_list[gain] = Symbol(gain)
        if plus_term not in self.var_list:
            self.var_list[plus_term] = Symbol(plus_term)
        if minus_term not in self.var_list:
            self.var_list[minus_term] = Symbol(minus_term)
        if ref_plus not in self.var_list:
            self.var_list[ref_plus] = Symbol(ref_plus)
        if ref_minus not in self.var_list:
            self.var_list[ref_minus] = Symbol(ref_minus)

        if plus_term not in self.nodes:
            self.nodes[plus_term] = []
        if minus_term not in self.nodes:
            self.nodes[minus_term] = []

        self.nodes[plus_term] += [StoredOp(ref_plus, ref_minus, gain, top=True, neg=False)]
        self.nodes[minus_term] += [StoredOp(ref_plus, ref_minus, gain, top=True, neg=True)]

    def _add_isource(self, plus_term, minus_term):
        if plus_term not in self.var_list:
            self.var_list[plus_term] = Symbol(plus_term)
        if minus_term not in self.var_list:
            self.var_list[minus_term] = Symbol(minus_term)

        if plus_term not in self.nodes:
            self.nodes[plus_term] = []
        if minus_term not in self.nodes:
            self.nodes[minus_term] = []

        if "Isource" not in self.var_list.keys():
            self.var_list["Isource"] = Symbol("Isource")

        self.nodes[plus_term] += [StoredOp(self.var_list["Isource"], neg=False, isource=True)]
        self.nodes[minus_term] += [StoredOp(self.var_list["Isource"], neg=True, isource=True)]

    def _prep_solve(self):
        remove_list = []
        for k in self.nodes.keys():
            node = self.nodes[k]
            if len(node) == 1 or k == "Gnd":
                remove_list += [k]
                self.var_list[k] = 0
        for k in remove_list:
            self.nodes.pop(k)
            self.var_list.pop(k)

        for k in self.nodes.keys():
            self.nodes[k] = sum([x.construct(self.var_list) for x in self.nodes[k]])

    def _save_state(self):
        self._bck_nodes = copy.deepcopy(self.nodes)
        self._bck_var_list = copy.deepcopy(self.var_list)

    def _finish_solve(self):
        self.nodes = self._bck_nodes
        self._bck_nodes = None
        self.var_list = self._bck_var_list
        self._bck_var_list = None

    def get_tf(self, in_name, out_name):
        in_name = Circ.format_name(in_name)
        out_name = Circ.format_name(out_name)
        if not in_name in self.var_list:
            raise ValueError("Input Node \'" + str(in_name) + "\' Does Not Exist")
        if not out_name in self.var_list:
            raise ValueError("Output Node" + str(out_name) + " Does Not Exist")

        self._prep_solve()

        solve_nodes = self._get_solve_path(in_name, out_name)
        replace_dict = dict()
        for i in range(len(solve_path) - 1):
            self.nodes[solve_path[i]] = self.nodes[solve_path[i]].xreplace(replace_dict)
            replace_dict[self.var_list[solve_path[i + 1]]] = solve(self.nodes[solve_path[i]], self.var_list[solve_path[i + 1]])[0]

        if mode == "i":
            tf = simplify((self.var_list[out_name] / self.nodes[in_name]).xreplace(replace_dict))
        elif mode == "v":
            tf = simplify(replace_dict[self.var_list[in_name]] / self.var_list[out_name])

        self._finish_solve()
        return tf

    def get_imp(self, vout_plus, vout_minus="Gnd", iin_plus=None, iin_minus="Gnd"):
        if iin_plus == None:
            iin_plus = vout_plus
        iin_plus = Circ.format_name(iin_plus)
        iin_minus = Circ.format_name(iin_minus)
        vout_plus = Circ.format_name(vout_plus)
        vout_minus = Circ.format_name(vout_minus)
        nodes_list = list(self.nodes.keys())

        if vout_plus not in nodes_list:
            vout_plus = "Gnd"
        if vout_minus not in nodes_list:
            vout_minus = "Gnd"
        if iin_plus not in nodes_list:
            iin_plus = "Gnd"
        if iin_minus not in nodes_list:
            iin_minus = "Gnd"

        if not vout_plus in self.var_list:
            raise ValueError("Input Node \'" + str(vout_plus) + "\' Does Not Exist")
        if not vout_minus in self.var_list:
            raise ValueError("Input Node \'" + str(vout_minus) + "\' Does Not Exist")
        if not iin_plus in self.var_list:
            raise ValueError("Input Node \'" + str(iin_plus) + "\' Does Not Exist")
        if not iin_minus in self.var_list:
            raise ValueError("Input Node \'" + str(iin_minus) + "\' Does Not Exist")

        self._save_state()
        self._add_isource(iin_plus, iin_minus)
        self._prep_solve()

        solve_nodes = list(self.nodes.keys())
        # solve_nodes = self._get_solve_nodes([n for n in [iin_plus, iin_minus, vout_plus, vout_minus] if n != "Gnd"])
        solve_nodes = self._get_solve_nodes_naive([n for n in [iin_plus, iin_minus, vout_plus, vout_minus] if n != "Gnd"])
        all_eq_list = []
        remainder_list = []
        for cur_node_name in solve_nodes:
            cur_eq_list = []
            cur_expr = self.nodes[cur_node_name].expand()
            cur_sum = 0
            for cur_coeff_name in solve_nodes:
                cur_coeff = cur_expr.coeff(self.var_list[cur_coeff_name])
                cur_eq_list += [cur_coeff]
                cur_sum += cur_coeff * self.var_list[cur_coeff_name]
            remainder = simplify(cur_expr - cur_sum)
            remainder_list += [remainder]
            all_eq_list += [cur_eq_list]

        a_matrix = Matrix(all_eq_list)
        b_matrix = Matrix(remainder_list)
        soln = a_matrix.LUsolve(b_matrix)
        vplus = soln[solve_nodes.index(vout_plus)] if vout_plus != "Gnd" else 0
        vminus = soln[solve_nodes.index(vout_minus)] if vout_minus != "Gnd" else 0
        tf = simplify((vplus - vminus) / self.var_list["Isource"])
        self._finish_solve()
        return tf

    def _get_solve_path(self, start, end):
        cur_graph = self._build_graph()
        cur_node = cur_graph[end]
        cur_node_name = end
        traversed_nodes = [end]
        while cur_node_name != start:
            unchanged = True
            for v in cur_node:
                if v not in traversed_nodes:
                    if start == v and len(traversed_nodes) < len(cur_node):
                        pass
                    else:
                        unchanged = False
                        break
            if unchanged:
                v = traversed_nodes[traversed_nodes.index(cur_node_name) - 1]
                cur_node = cur_graph[v]
                cur_node_name = v
            else:
                traversed_nodes += [v]
                cur_node = cur_graph[v]
                cur_node_name = v
        return traversed_nodes

    def _get_solve_nodes_naive(self, nodes_list):
        all_symbols = set(chain(*[self._get_expr_nodes(n) for n in nodes_list]))
        nodes_list = set(nodes_list)
        while all_symbols != nodes_list:
            nodes_list = all_symbols
            all_symbols = set(chain(*[self._get_expr_nodes(n) for n in nodes_list]))
        return list(nodes_list)

    def _get_solve_nodes(self, nodes_list):
        # all_symbols = set(chain([self._get_expr_nodes(self.nodes[n]) for n in nodes_list]**))
        num_eq = len(nodes_list)
        num_symbols = len(all_symbols)
        solve_nodes = nodes_list
        cur_node_name = nodes_list[0]
        cur_node = self.nodes[cur_node_name]
        while num_eq < num_symbols:
            for v in self._get_expr_nodes(cur_node):
                if v not in solve_nodes:
                    overlap, _, extra_newCirc._num_overlap(self._get_expr_nodes[cur_node], all_symbols)
                    if overlap < 2:
                        pass
                    else:
                        unchanged = False
                        break
            if unchanged:
                # v = traversed_nodes[traversed_nodes.index(cur_node_name) + 1]
                # cur_node = cur_graph[v]
                # cur_node_name = v
                pass
            else:
                # solve_nodes +=
                pass
        return solve_nodes

    def _get_expr_nodes(self, node_name, ret_all=False):
        expr = self.nodes[node_name]
        nodes = [str(x) for x in expr.free_symbols if "V" == str(x)[0] or "v" == str(x)[0]]
        if node_name in nodes or ret_all:
            return nodes
        return []

    def _num_overlap(lst_a, lst_b, ret_all=False):
        count = 0
        for elem in lst_a:
            if elem in lst_b:
                count += 1
        return count, len(lst_a) - count, len(lst_b) - count

    def _get_imp_path(self, in_name):
        cur_graph = self._build_graph()
        cur_node = cur_graph[in_name]
        cur_node_name = in_name
        traversed_nodes = [in_name]
        while True:
            unchanged = True
            for v in cur_node:
                if v not in traversed_nodes:
                    unchanged = False
                    break
            if unchanged:
                return traversed_nodes
            traversed_nodes += [v]
            cur_node = cur_graph[v]
            cur_node_name = v


    def _build_graph(self):
        leaf_list = dict()
        for k in self.nodes.keys():
            leaf_list[k] = []
            for v in self.var_list.keys():
                if self.nodes[k].has(v) and v != k and "V" in v:
                    leaf_list[k] += [v]
        return leaf_list

    def print_contents(self):
        if type(self.nodes[list(self.nodes.keys())[0]]) == list:
            mod = True
            self._save_state()
            self._prep_solve()
        else:
            mod = False
        print("NODES INFO\n")
        for k in self.nodes.keys():
            print("Node", k)
            pprint.pprint(self.nodes[k])
            print()

        print()
        print("Variables")
        pprint.pprint(list(self.var_list.keys()))
        if mod:
            self._finish_solve()

    def format_name(s):
        c = s[0]
        rest = s[1:]
        if not c.isalpha():
            raise ValueError("Entered name " + s + " is not valid. First character is not part of the alphabet")
        return c.upper() + rest.lower()

class StoredOp():

    def __init__(self, plus="", minus="", factor=1, top=True, add_jw=False, neg=False, isource=False):
        self.plus = plus
        self.minus = minus
        self.factor = factor
        self.top = top
        self.add_jw = add_jw
        self.neg = neg
        self.isource = isource

    def construct(self, var_list):
        if self.isource:
            if self.neg:
                return -var_list["Isource"]
            else:
                return var_list["Isource"]

        plus = var_list[self.plus] if self.plus != "Gnd" else 0
        minus = var_list[self.minus] if self.minus != "Gnd" else 0
        diff = plus - minus
        if self.neg:
            diff *= -1

        custom_expr = var_list[self.factor]
        if self.add_jw:
            custom_expr *= Circ.jw
        if self.top:
            new_expr = diff * custom_expr
        else:
            new_expr = diff / custom_expr

        return new_expr


if __name__ == "__main__":
    # Rfb = Symbol("Rfb")
    # Gm = Symbol("Gm")
    # Cpd = Symbol("Cpd")
    # L = Symbol("L")
    # Ro = Symbol("Ro")
    # Vout = Symbol("Vout")
    # Vin = Symbol("Vin")
    # Vmid = Symbol("Vmid")
    # Gnd = Symbol("Gnd")
    # jw = Symbol("jw")

    circ = Circ()
    circ.add_res("rfb", "vin", "vout")
    circ.add_cond("gds", "vout", "gnd")
    circ.add_vccs("gm", "vout", "gnd", "gnd", "Vin")
    circ.add_cap("cpd", "vin", "gnd")

    # circ.print_contents()
    tf = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    nf = circ.get_imp(vout_plus="Vout", iin_plus="Vout")
    print("tf", tf)
    print("nf", nf)
    print('nf/tf', nf / tf)

    # print("Test 2")
    # circ = Circ()
    # circ.add_res("rfb", "vmid", "vout")
    # circ.add_res("ro", "gnd", "vout")
    # circ.add_vccs("gm", "vout", "gnd", "vmid", "gnd")
    # circ.add_cap("cpd", "vin", "gnd")
    # circ.add_ind("l", "vin", "vmid")

    # circ.print_contents()
    # imp = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    # print(imp)
    # print()

    # print("Test 3")
    # circ = Circ()
    # circ.add_cap("cpd", "vin", "gnd")
    # circ.add_res("rfb", "vout", "vin")
    # circ.add_vccs("gm", "gnd", "vmid", "vin", "gnd")
    # circ.add_vccs("gds", "gnd", "Vmid", "vmid", "gnd")
    # circ.add_vccs("A", "gnd", "Vout", "Vmid", "gnd")
    # circ.add_vccs("gout", "Vout", "gnd", "gnd", "vout")
    # tf = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    # nf = circ.get_imp(vout_plus="Vout", iin_plus="Vmid")
    # print(imp)

    print("Test 4")
    circ = Circ()
    circ.add_cap("cpd", "vin", "gnd")
    circ.add_res("rfb", "vout", "vin")
    circ.add_vccs("gm_inv", "gnd", "vmid", "vin", "gnd")
    circ.add_vccs("gds_inv", "gnd", "Vmid", "vmid", "gnd")
    circ.add_vccs("gm_boost", "gnd", "Vout", "Vmid", "gnd")
    circ.add_vccs("gds_boost", "Vout", "gnd", "gnd", "vout")
    tf = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    nf = circ.get_imp(vout_plus="Vout", iin_plus="Vmid")

