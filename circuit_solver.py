from sympy import Symbol, solve, simplify
import copy
import pprint

class Circ():
    """
    All nodes and values are not case Sensitive (defaults to Capitol followed by lower case)
    Gnd, jw, kbt: reserved names

    Voltages start with V
    Inductors start with L
    Resistors start with R
    Capacitors start with C

    """

    jw = Symbol("jw")
    kbt = Symbol("kbt")

    def __init__(self):
        self.nodes = dict()
        self.var_list = dict()


    def add_res(self, res, plus_term, minus_term):
        res = format_string(res)
        plus_term = format_string(plus_term)
        minus_term = format_string(minus_term)

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

        # self.nodes[a] += [(self.var_list[a] - self.var_list[b]) / self.var_list[res]]
        # self.nodes[b] += [-(self.var_list[a] - self.var_list[b]) / self.var_list[res]]
        self.nodes[plus_term] += [StoredOp(plus_term, minus_term, res, top=False, add_jw=False, neg=True)]
        self.nodes[minus_term] += [StoredOp(plus_term, minus_term, res, top=False, add_jw=False, neg=False)]


    def add_cap(self, cap, plus_term, minus_term):
        cap = format_string(cap)
        plus_term = format_string(plus_term)
        minus_term = format_string(minus_term)

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

        # self.nodes[a] += [(self.var_list[a] - self.var_list[b]) * self.var_list[cap] * Circ.jw]
        # self.nodes[b] += [-(self.var_list[a] - self.var_list[b]) * self.var_list[cap] * Circ.jw]
        self.nodes[plus_term] += [StoredOp(plus_term, minus_term, cap, top=True, add_jw=True, neg=True)]
        self.nodes[minus_term] += [StoredOp(plus_term, minus_term, cap, top=True, add_jw=True, neg=False)]


    def add_ind(self, ind, plus_term, minus_term):
        ind = format_string(ind)
        plus_term = format_string(plus_term)
        minus_term = format_string(minus_term)

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

        # self.nodes[a] += [(self.var_list[a] - self.var_list[b]) / (self.var_list[ind] * Circ.jw)]
        # self.nodes[b] += [-(self.var_list[a] - self.var_list[b]) / (self.var_list[ind] * Circ.jw)]
        self.nodes[plus_term] += [StoredOp(plus_term, minus_term, ind, top=False, add_jw=True, neg=True)]
        self.nodes[minus_term] += [StoredOp(plus_term, minus_term, ind, top=False, add_jw=True, neg=False)]


    def add_vcvs(self, gain, plus_term, minus_term, ref_plus, ref_minus):
        gain = format_string(gain)
        plus_term = format_string(plus_term)
        minus_term = format_string(minus_term)
        ref_plus = format_string(ref_plus)
        ref_minus = format_string(ref_minus)

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

        # self.nodes[a] += [self.var_list[gain] * (self.var_list[ref_plus] - self.var_list[ref_minus])]
        # self.nodes[b] += [-self.var_list[gain] * (self.var_list[ref_plus] - self.var_list[ref_minus])]
        self.nodes[plus_term] += [StoredOp(ref_plus, ref_minus, gain, top=True, neg=False)]
        self.nodes[minus_term] += [StoredOp(ref_plus, ref_minus, gain, top=True, neg=True)]


    def _prep_solve(self):
        # for k in self.nodes.keys():
        #     self.nodes[k] = [x.construct(self.var_list) for x in self.nodes[k]]
        # pprint.pprint(self.nodes)
        # return
        self._bck_nodes = copy.deepcopy(self.nodes)
        self._bck_var_list = copy.deepcopy(self.var_list)
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


    def _finish_solve(self):
        self.nodes = self._bck_nodes
        self._bck_nodes = None
        self.var_list = self._bck_var_list
        self._bck_var_list = None


    def get_tf(self, in_name, out_name, mode="v"):
        in_name = format_string(in_name)
        out_name = format_string(out_name)
        if not in_name in self.var_list:
            raise ValueError("Input Node \'" + str(in_name) + "\' Does Not Exist")
        if not out_name in self.var_list:
            raise ValueError("Output Node" + str(out_name) + " Does Not Exist")
        orig_mode = mode
        mode = mode.lower()
        if mode not in ["v", "i"]:
            raise ValueError("Invalid Mode \'" + str(orig_mode) + "\'")

        self._prep_solve()

        solve_path = self._get_solve_path(in_name, out_name)
        print(solve_path)
        replace_dict = dict()
        # ref_node = solve_path[0]
        for i in range(len(solve_path) - 1):
            self.nodes[solve_path[i]] = self.nodes[solve_path[i]].xreplace(replace_dict)
            replace_dict[self.var_list[solve_path[i + 1]]] = solve(self.nodes[solve_path[i]], self.var_list[solve_path[i + 1]])[0]

        if mode == "i":
            tf = simplify((self.var_list[out_name] / self.nodes[in_name]).xreplace(replace_dict))
        elif mode == "v":
            tf = simplify(replace_dict[self.var_list[in_name]] / self.var_list[out_name])

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
                    v = missing_node
            traversed_nodes += [v]
            cur_node = cur_graph[v]
            cur_node_name = v
        return traversed_nodes

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


    def get_imp(self, in_name):
        in_name = format_string(in_name)
        if not in_name in self.var_list:
            raise ValueError("Input Node \'" + str(in_name) + "\' Does Not Exist")

        self._prep_solve()

        solve_path = self._get_imp_path(in_name)[::-1]
        replace_dict = dict()
        print(solve_path)

        for i in range(len(solve_path) - 1):
            x = solve(self.nodes[solve_path[i]].xreplace(replace_dict), self.var_list[solve_path[i]])[0]
            print(solve_path[i], x)
            replace_dict[self.var_list[solve_path[i]]] = solve(self.nodes[solve_path[i]].xreplace(replace_dict), self.var_list[solve_path[i]])[0]
            for k in replace_dict.keys():
                replace_dict[k] = replace_dict[k].xreplace(replace_dict)


        tf = simplify(self.var_list[in_name] / (self.nodes[in_name].xreplace(replace_dict)))
        self._finish_solve()
        return tf

    def print_contents(self):
        # DONT CALL FOR NOW
        if type(self.nodes[list(self.nodes.keys())[0]]) == list:
            mod = True
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


class StoredOp():

    def __init__(self, plus, minus, factor, top=True, add_jw=False, neg=False):
        self.plus = plus
        self.minus = minus
        self.factor = factor
        self.top = top
        self.add_jw = add_jw
        self.neg = neg


    def construct(self, var_list):
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


def format_string(s):
    c = s[0]
    rest = s[1:]
    if not c.isalpha():
        raise ValueError("Entered name " + s + " is not valid. First character is not part of the alphabet")
    return c.upper() + rest.lower()

def eval_func(tf_expr, val_dict):
    pass

l = locals()

def make_var(circ):
    save_dict = dict()
    l["_var_list"] = list(circ.var_list.keys())
    for v in circ.var_list.keys():
        if v in l:
            save_dict[v] = l[v]
        l[v] = Symbol(v)
    l["_save_dict"] = save_dict

def clear_var():
    pass
if __name__ == "__main__":
    Rfb = Symbol("Rfb")
    Gm = Symbol("Gm")
    Cpd = Symbol("Cpd")
    L = Symbol("L")
    Ro = Symbol("Ro")
    Vout = Symbol("Vout")
    Vin = Symbol("Vin")
    Vmid = Symbol("Vmid")
    Gnd = Symbol("Gnd")
    jw = Symbol("jw")

    circ = Circ()
    circ.add_res("rfb", "vin", "vout")
    circ.add_res("ro", "gnd", "vout")
    circ.add_vcvs("gm", "vout", "gnd", "vin", "gnd")
    circ.add_cap("cpd", "vin", "gnd")

    # circ.print_contents()
    tf = circ.get_tf("Vin", "Vout", "i")
    print(tf)
    print()

    print("Test 2")
    circ = Circ()
    circ.add_res("rfb", "vmid", "vout")
    circ.add_res("ro", "gnd", "vout")
    circ.add_vcvs("gm", "vout", "gnd", "vmid", "gnd")
    circ.add_cap("cpd", "vin", "gnd")
    circ.add_ind("l", "vin", "vmid")

    # circ.print_contents()
    # tf = circ.get_tf("Vin", "Vout", "i")
    print(tf)
    print()
    imp = circ.get_imp("Vout")
    print(imp)
    print()


    print("Test 3")
    circ = Circ()
    circ.add_cap("cpd", "vin", "gnd")
    circ.add_res("rfb", "vout", "vin")
    circ.add_vcvs("gm", "gnd", "vmid", "vin", "gnd")
    circ.add_vcvs("gds", "gnd", "Vmid", "vmid", "gnd")
    circ.add_vcvs("A", "gnd", "Vout", "Vmid", "gnd")
    circ.add_vcvs("gout", "Vout", "gnd", "gnd", "vout")
    tf = circ.get_tf("Vin", "Vout", "i")
    nf = circ.get_tf("Vmid", "Vout", "i")
    print(tf)
    print(imp)


