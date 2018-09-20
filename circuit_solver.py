from sympy import Symbol, solve, simplify, expand, Matrix
import sympy
import copy
import pprint
import numbers
from itertools import chain

class Circ():
    """Main Circuit Solver Class

    All nodes and values are not case Sensitive (defaults to Capitol followed by lower case)
    Gnd, jw, kbt, Isource: reserved names

    Voltages should start with V
    Inductors should start with L
    Resistors should start with R
    Capacitors should start with C

    """
    jw = Symbol("jw")
    w = Symbol("w")
    j = Symbol("j")
    kbt = Symbol("kbt")

    def __init__(self):
        self.nodes = dict()
        self.var_list = dict()
        self.exclude_list = set()

    def clear(self):
        self.nodes = dict()
        self.var_list = dict()
        self.exclude_list = set()

    def add_res(self, res, plus_term, minus_term):
        """Add a resistor
        """
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

        self.nodes[plus_term] += [_StoredOp(plus_term, minus_term, res, top=False, add_jw=False, neg=True)]
        self.nodes[minus_term] += [_StoredOp(plus_term, minus_term, res, top=False, add_jw=False, neg=False)]

    def add_conductance(self, cond, plus_term, minus_term):
        """Add a conductance
        """
        self.add_vccs(cond, plus_term, minus_term, minus_term, plus_term)

    def add_cap(self, cap, plus_term, minus_term):
        """Add a capacitor
        """
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

        self.nodes[plus_term] += [_StoredOp(plus_term, minus_term, cap, top=True, add_jw=True, neg=True)]
        self.nodes[minus_term] += [_StoredOp(plus_term, minus_term, cap, top=True, add_jw=True, neg=False)]

    def add_ind(self, ind, plus_term, minus_term):
        """Add an inductor
        """
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

        self.nodes[plus_term] += [_StoredOp(plus_term, minus_term, ind, top=False, add_jw=True, neg=True)]
        self.nodes[minus_term] += [_StoredOp(plus_term, minus_term, ind, top=False, add_jw=True, neg=False)]

    def add_vccs(self, gain, plus_term, minus_term, ref_plus, ref_minus):
        """Add a voltage controlled current source (useful for simulating transistors)
        """
        gain = Circ.format_name(gain)
        plus_term = Circ.format_name(plus_term)
        minus_term = Circ.format_name(minus_term)
        ref_plus = Circ.format_name(ref_plus)
        ref_minus = Circ.format_name(ref_minus)
        if ref_plus != "Gnd":
            self.exclude_list.add(ref_plus)
        if ref_minus != "Gnd":
            self.exclude_list.add(ref_minus)

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
        if ref_plus not in self.nodes:
            self.nodes[ref_plus] = []
        if ref_minus not in self.nodes:
            self.nodes[ref_minus] = []

        self.nodes[plus_term] += [_StoredOp(ref_plus, ref_minus, gain, top=True, neg=False)]
        self.nodes[minus_term] += [_StoredOp(ref_plus, ref_minus, gain, top=True, neg=True)]

    def add_transistor(self, vd='Vd', vg='Vg', vs='Gnd', gm=None, gds=None, ro=None, cl=None):
        """Add a small signal model of a transistor formed by 'gm', 'gds' (or 'ro') and a potentially a load capacitor 'cl'
        """
        if gm == None:
            raise ValueError("Missing gm value")
        if gds == None and ro == None:
            raise ValueError("Missing ro and gds value")
        if gds != None:
            self.add_conductance(gds, vd, vs)
        else:
            self.add_res(ro, vd, vs)
        self.add_vccs(gm, vs, vd, vg, vs)
        if cl == None:
            return
        self.add_cap(cl, vd, "Gnd")

    def get_tf(self, vout_plus, vin_plus, vout_minus="Gnd", vin_minus="Gnd"):
        """Caluclate the Voltage Transfer Function
        """
        if vin_plus == None:
            vin_plus = vout_plus
        vin_plus = Circ.format_name(vin_plus)
        vin_minus = Circ.format_name(vin_minus)
        vout_plus = Circ.format_name(vout_plus)
        vout_minus = Circ.format_name(vout_minus)
        nodes_list = list(self.nodes.keys())

        if not vout_plus in self.var_list:
            raise ValueError("Input Node \'" + str(vout_plus) + "\' Does Not Exist")
        if not vout_minus in self.var_list:
            raise ValueError("Input Node \'" + str(vout_minus) + "\' Does Not Exist")
        if not vin_plus in self.var_list:
            raise ValueError("Input Node \'" + str(vin_plus) + "\' Does Not Exist")
        if not vin_minus in self.var_list:
            raise ValueError("Input Node \'" + str(vin_minus) + "\' Does Not Exist")

        self._save_state()
        self._prep_solve()
        nodes_list_solve = self.nodes
        self._finish_solve()

        self._save_state()
        self._prep_solve(exclude_list=self.exclude_list)

        in_nodes = [n for n in [vin_plus, vin_minus] if n != "Gnd"]
        out_nodes = [n for n in [vout_plus, vout_minus] if n != "Gnd"]
        solve_nodes = self._get_tf_nodes(in_nodes=in_nodes, out_nodes=out_nodes)
        solve_nodes = [n for n in solve_nodes if n in nodes_list_solve]
        soln = self._solve_matrix(solve_nodes)

        vout_plus_val = soln[solve_nodes.index(vout_plus)] if vout_plus != "Gnd" else 0
        vout_minus_val = soln[solve_nodes.index(vout_minus)] if vout_minus != "Gnd" else 0
        vin_plus_val = self.var_list[vin_plus] if vin_plus != "Gnd" else 0
        vin_minus_val = self.var_list[vin_minus] if vin_minus != "Gnd" else 0
        tf = simplify((vout_plus_val - vout_minus_val) / (vin_plus_val - vin_minus_val))
        self._finish_solve()
        return tf

    def get_imp(self, vout_plus, vout_minus="Gnd", iin_plus=None, iin_minus="Gnd"):
        """Calucate the Impedance Transfer Function
        """
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

        if iin_plus not in self.nodes and iin_minus not in self.nodes:
            # raise ValueError("Impossible to calculate impedance transfer function. Output voltage is independent of input current. (Most likely infinite impedence)")
            print("Impossible to calculate transfer function. Output voltage is independent of input current. Evaluating as infinite impedence)")
            self._finish_solve()
            return float('inf')

        solve_nodes = self._get_imp_nodes([n for n in [iin_plus, iin_minus, vout_plus, vout_minus] if n != "Gnd"])
        soln = self._solve_matrix(solve_nodes)

        if (vout_plus not in solve_nodes and vout_plus != "Gnd") or (vout_minus not in solve_nodes and vout_minus != "Gnd"):
            # raise ValueError("Impossible to calculate transfer function. Output voltage is independent of input current. (Most likely infinite impedence)")
            print("Impossible to calculate impedance transfer function. Output voltage is independent of input current. Evaluating as infinite impedence)")
            self._finish_solve()
            return float('inf')

        vplus = soln[solve_nodes.index(vout_plus)] if vout_plus != "Gnd" else 0
        vminus = soln[solve_nodes.index(vout_minus)] if vout_minus != "Gnd" else 0
        tf = simplify((vplus - vminus) / self.var_list["Isource"])
        tf = simplify(expand(tf))
        self._finish_solve()
        return tf

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

        self.nodes[plus_term] += [_StoredOp(self.var_list["Isource"], neg=False, isource=True)]
        self.nodes[minus_term] += [_StoredOp(self.var_list["Isource"], neg=True, isource=True)]

    def _prep_solve(self, exclude_list=set()):
        remove_list = []
        for k in self.nodes.keys():
            node = self.nodes[k]
            if (len(node) == 1 or k == "Gnd") and k not in exclude_list:
                remove_list += [k]
                self.var_list[k] = 0
            elif (len(node) == 0 or k == "Gnd") and k not in exclude_list:
                remove_list += [k]
        for k in remove_list:
            self.nodes.pop(k)

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

    def _get_tf_nodes(self, in_nodes, out_nodes):
        """
        in: Vin
        out: Vout

        Vout: (Vout, Vin)
        """
        nodes_list = out_nodes
        all_symbols = set(chain(*[self._get_expr_nodes(n) for n in nodes_list]))
        while all_symbols != set(list(nodes_list) + in_nodes):
            nodes_list = all_symbols
            all_symbols = set(chain(*[self._get_expr_nodes(n) for n in nodes_list if n in self.nodes and not isinstance(self.nodes[n], numbers.Number)]))
        return [n for n in list(nodes_list) if n in self.nodes]

    def _get_imp_nodes(self, nodes_list):
        all_symbols = set(chain(*[self._get_expr_nodes(n) for n in nodes_list]))
        while all_symbols != nodes_list:
            nodes_list = all_symbols
            all_symbols = set(chain(*[self._get_expr_nodes(n) for n in nodes_list]))
        return list(nodes_list)

    def _solve_matrix(self, solve_nodes):
        all_eq_list = []
        remainder_list = []
        for cur_node_name in solve_nodes:
            cur_eq_list = []
            cur_expr = self.nodes[cur_node_name].expand() if not isinstance(self.nodes[cur_node_name], numbers.Number) else self.nodes[cur_node_name]
            cur_sum = 0
            for cur_coeff_name in solve_nodes:
                cur_coeff = cur_expr.coeff(self.var_list[cur_coeff_name])
                cur_eq_list += [cur_coeff]
                cur_sum += cur_coeff * self.var_list[cur_coeff_name]
            remainder = simplify(cur_sum - cur_expr)
            remainder_list += [remainder]
            all_eq_list += [cur_eq_list]

        a_matrix = Matrix(all_eq_list)
        b_matrix = Matrix(remainder_list)
        soln = a_matrix.LUsolve(b_matrix)
        return soln

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

    def print_contents(self):
        """Print information about the modeled circuit
        """
        if type(self.nodes[list(self.nodes.keys())[0]]) == list:
            mod = True
            self._save_state()
            self._prep_solve(exclude_list=self.exclude_list)
        else:
            mod = False
        print("NODES INFO")
        for k in self.nodes.keys():
            print("Node", k)
            pprint.pprint(self.nodes[k])
            print()

        print("Variables")
        pprint.pprint(list(self.var_list.keys()))
        if mod:
            self._finish_solve()

    @property
    def exec_cmd(self):
        """Make a command that, when called, sets all the used variables to be symbols
        This way it is much easier to play around with the transfer functions
        """
        exec_str = ""
        for var in set(list(self.var_list.keys()) + ['j', 'w', 'jw', 'kbt', 'Gnd']):
            exec_str += var + "=Symbol(\'" + var + "\')\n"
        return exec_str

    def format_name(s):
        """Formats a string such that the first letter is capitolized and the rest lower cased
        """
        c = s[0]
        rest = s[1:]
        if not c.isalpha():
            raise ValueError("Entered name " + s + " is not valid. First character is not part of the alphabet")
        return c.upper() + rest.lower()

class _StoredOp():
    """Helper Class for representing the circuit
    """
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
    pass

    """
    some random test cases and old code
    """

    # print("Test 1")
    # circ = Circ()
    # circ.add_res("rfb", "vin", "vout")
    # circ.add_conductance("gds", "vout", "gnd")
    # circ.add_vccs("gm", "vout", "gnd", "gnd", "Vin")
    # circ.add_cap("cpd", "vin", "gnd")

    # circ.print_contents()
    # tf1 = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    # nf1 = circ.get_imp(vout_plus="Vout", iin_plus="Vout")
    # print("tf", tf1)
    # print("nf", nf1)
    # print('nf/tf', nf1 / tf1)

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

    # print("Test 4")
    # circ = Circ()
    # circ.add_cap("cpd", "vin", "gnd")
    # circ.add_res("rfb", "vout", "vin")
    # circ.add_vccs("gm_inv", "gnd", "vmid", "vin", "gnd")
    # circ.add_vccs("gds_inv", "gnd", "Vmid", "vmid", "gnd")
    # circ.add_vccs("gm_boost", "gnd", "Vout", "Vmid", "gnd")
    # circ.add_vccs("gds_boost", "Vout", "gnd", "gnd", "vout")
    # tf = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    # nf = circ.get_imp(vout_plus="Vout", iin_plus="Vmid")

    # print("Test 5")
    # circ = Circ()
    # circ.add_vccs("g1", "V2", "gnd", "V1", "gnd")
    # circ.add_vccs("g2", "V2", "gnd", "V3", "gnd")
    # circ.add_res("r1", "v1", "v3")
    # circ.add_res("r2", "v1", "gnd")
    # circ.add_res("r3", "gnd", "v3")

    # print("Test 6")
    # circ = Circ()
    # circ.add_cap("cpd", "vin", "gnd")
    # circ.add_res("rfb", "vout", "vin")
    # circ.add_vccs("gm", "gnd", "vout", "vin", "gnd")
    # circ.add_conductance("gds", "vout", "gnd")
    # circ.add_cap("cload", "vout", "gnd")
    # tf_orig = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    # circ.clear()

    # circ.add_cap("cpd", "vin", "gnd")
    # circ.add_res("rfb", "vout", "vin")
    # circ.add_vccs("gm", "gnd", "vmid", "vin", "gnd")
    # circ.add_conductance("gds", "vmid", "gnd")
    # circ.add_res("rwire", "vout", "vmid")
    # circ.add_cap("cload", "vout", "gnd")
    # tf_res_top = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    # circ.clear()

    # circ = Circ()
    # circ.add_cap("cpd", "vin", "gnd")
    # circ.add_res("rfb", "vout", "vmid")
    # circ.add_vccs("gm", "gnd", "vout", "vmid", "gnd")
    # circ.add_conductance("gds", "vout", "gnd")
    # circ.add_res("rwire", "vmid", "vin")
    # circ.add_cap("cload", "vout", "gnd")
    # tf_res_bot = circ.get_imp(vout_plus="Vout", iin_plus="Vin")

    # print("Test 7")
    # circ = Circ()
    # circ.add_vccs("Gm1", "Gnd", "Vout", "Vin", "Gnd")
    # circ.add_vccs("Gm2", "Gnd", "Vin", "Vout", "Gnd")
    # circ.add_conductance("Gds1", "Gnd", "Vout")
    # circ.add_conductance("Gds2", "Gnd", "Vin")
    # circ.add_res("Rfb", "Vin", "Vout")
    # circ.add_cap("Cpd", "Vin", "Gnd")
    # imp = circ.get_imp(vout_plus="Vout", iin_plus="Vin")
    # tf = circ.get_tf(vout_plus="Vout", vin_plus="Vin")

    # print("Test 8")
    # circ = Circ()
    # circ.add_vccs("Gm1", "Gnd", "Vmid", "Vin", "Gnd")
    # circ.add_vccs("Gm2", "Gnd", "Vout", "Vmid", "Gnd")
    # circ.add_conductance("Gds1", "Gnd", "Vout")
    # circ.add_conductance("Gds2", "Gnd", "Vmid")
    # tf = circ.get_tf(vout_plus="Vout", vin_plus="Vin")
    # tf_mid = circ.get_tf(vout_plus="Vmid", vin_plus="Vin")

    # print("Test 9")
    # circ = Circ()
    # circ.add_transistor(vd='Vd', vg='Vg', vs='Gnd', gm='gm', gds='gds', cl='cl')
    # trans_imp = circ.get_imp(vout_plus="Vd", iin_plus="Vg")
    # trans_tf = circ.get_tf(vout_plus="Vd", vin_plus="Vg")

    # print("Test 10")
    # circ = Circ()
    # circ.add_res("Rj", "Vin", "Gnd")
    # circ.add_cap("Cj", "Vin", "Gnd")
    # circ.add_res("Rs", "Vin", "V1")
    # circ.add_cap("Cp", "V1", "Gnd")
    # circ.add_ind("Ls", "V1", "V2")
    # circ.add_vccs("Gm", "Gnd", "Vout", "V2", "Gnd")
    # circ.add_conductance("Gds", "Vout", "Gnd")
    # circ.add_res("Rfb", "Vout", "V2")
    # imp = circ.get_imp(vout_plus="Vout", iin_plus="Vin")

    # print("Test 11")
    # circ = Circ()
    # circ.add_res("Rj", "Vin", "Gnd")
    # circ.add_cap("Cj", "Vin", "Gnd")
    # circ.add_res("Rs", "Vin", "V1")
    # circ.add_cap("Cp", "V1", "Gnd")
    # circ.add_res("Zin", "V1", "Gnd")
    # circ.add_vccs("Gm", "Gnd", "Vout", "V1", "Gnd")
    # circ.add_conductance("Gds", "Vout", "Gnd")
    # circ.add_res("Rfb", "Vout", "V1")
    # imp = circ.get_imp(vout_plus="V1", iin_plus="Vin")
    # imp = circ.get_imp(vout_plus="Vout", iin_plus="Vin")

    # Cj = 250e-15
    # Cp = 20e-15
    # Rj = 1e6
    # Rs = 5
    # jw = Symbol("jw")


    # print("Test 12")
    # circ = Circ()
    # circ.add_res("rfb", "vin", "vout")
    # circ.add_conductance("gds", "vout", "gnd")
    # circ.add_vccs("gm", "vout", "gnd", "gnd", "Vin")
    # circ.add_cap("cpd", "vin", "gnd")
    # circ.add_res("rpd", "vin", "gnd")

    # tf_with = circ.get_imp(vout_plus="Vout", iin_plus="Vin")

    # circ.clear()
    # circ.add_res("rfb", "vin", "vout")
    # circ.add_conductance("gds", "vout", "gnd")
    # circ.add_vccs("gm", "vout", "gnd", "gnd", "Vin")
    # circ.add_cap("cpd", "vin", "gnd")

    # tf_wout = circ.get_imp(vout_plus="Vout", iin_plus="Vin")

