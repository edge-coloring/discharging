import argparse 
import os 
import re

INF = 10000

def deg2int(deg_str: str):
    if deg_str[-1] == '+' or deg_str[-1] == '-':
        return int(deg_str[:-1])
    else:
        return int(deg_str)
    
#
# 多重辺ありのときは
#   N E N_multiedge deg0 deg1 .. deg{N-1} u0 v0 u1 v1 .. u{E-1} v{E-1} u'0 v'0 u'1 v'1 .. u'{N_multiedge-1} v'{N_multiedge-1}
# 多重辺なしのときは
#   N E deg0 deg1 .. deg{N-1} u0 v0 u1 v1 .. u{E-1} v{E-1}
# の形の文字列を cartwheel の情報にして返す。
# 
def cartWheel(cartwheel_str: str, is_parallel: bool):
    if is_parallel:
        info_list = cartwheel_str.split(' ')
        N = int(info_list[0])
        E = int(info_list[1])
        N_multiedge = int(info_list[2])
        nodes = [i for i in range(N)]
        edges = []
        degrees = [0 for i in range(N)]
        show = [True for i in range(N)]
        multiedges = []
        for i in range(N):
            degi = info_list[3 + i]
            if degi == "?":
                show[i] = False
            elif degi == "9+":
                show[i] = False
            else:
                degrees[i] = deg2int(degi)
        for i in range(E):
            ui = int(info_list[3 + N + i * 2])
            vi = int(info_list[3 + N + i * 2 + 1])
            edges.append((ui, vi))
        for i in range(N_multiedge):
            ui = int(info_list[3 + N + 2 * E + i * 2])
            vi = int(info_list[3 + N + 2 * E + i * 2 + 1])
            multiedges.append((ui, vi))
        return N, nodes, edges, degrees, show, multiedges
    else:
        info_list = cartwheel_str.split(' ')
        N = int(info_list[0])
        E = int(info_list[1])
        nodes = [i for i in range(N)]
        edges = []
        degrees = [0 for i in range(N)]
        show = [True for i in range(N)]
        multiedges = []
        for i in range(N):
            degi = info_list[2 + i]
            if degi == "?":
                show[i] = False
            elif degi == "9+":
                show[i] = False
            else:
                degrees[i] = deg2int(degi)
        for i in range(E):
            ui = int(info_list[2 + N + i * 2])
            vi = int(info_list[2 + N + i * 2 + 1])
            edges.append((ui, vi))
        return N, nodes, edges, degrees, show, multiedges


#
# 頂点 v に対して N_v = {u | 辺 uv が存在する and show[u] = True } とする
# 頂点 v に対して M_v = {u | "多重" 辺 uv が存在する and show[u] = True } とする (M_v はサイズが 1 以下)
# + show[v] = true であるが、degrees[v] - |N_v| - |M_v| >= 4 以上の頂点 v 
# が存在しなくなるまで、そのような頂点 v を show[v] = false に変更する操作を続ける。
#
def markUnnecessaryVertices(N: int, edges: list, degrees: list, show: list, multiedges: list, is_remain: bool):
    incident = [set() for _ in range(N)]
    for e in edges:
        if show[e[0]]:
            incident[e[1]].add(e[0])
        if show[e[1]]:
            incident[e[0]].add(e[1])
    incident_multi = [set() for _ in range(N)]
    for e in multiedges:
        if show[e[0]]:
            incident_multi[e[1]].add(e[0])
        if show[e[1]]:
            incident_multi[e[0]].add(e[1])
    
    def dfs(v: int):
        deletion = False
        if is_remain:
            # is_remain == True のときは隣接している頂点が 1 つのときだけ消す。
            deletion = (len(incident[v]) == 1)
        else:
            # is_remain == False のときは 4 本以上の辺が外に出ているときに消す。
            deletion = (degrees[v] - len(incident[v]) - len(incident_multi[v]) >= 4)
        if show[v] == True and deletion:
            show[v] = False
        else:
            return
        for u in incident[v]:
            incident[u].remove(v)
        for u in incident_multi[v]:
            incident_multi[u].remove(v)
        for u in incident[v]:
            dfs(u)
        return
        
    for v in range(N):
        dfs(v)

    return


# diaognal[(u, v)] := {w : uvw が G 上で triangle になっている。} を計算する。
# configuration の場合は len(diagonal[e]) <= 2 を満たす。
def calcDiagonalVertices(N: int, edges: list):
    VtoV = [set() for _ in range(N)]
    diagonal = dict()
    for e in edges:
        VtoV[e[0]].add(e[1])
        VtoV[e[1]].add(e[0])
        diagonal[e] = []
    for e in diagonal.keys():
        for u in VtoV[e[0]]:
            if u in VtoV[e[1]]:
                diagonal[e].append(u)
        assert len(diagonal[e]) <= 2
    
    return diagonal


# カット点の列挙
def getCutVertices(N: int, VtoV: list):
    num = [-1 for _ in range(N)]
    low = [-1 for _ in range(N)]
    cutvertices = []
    def dfs(v: int, par: int, order: int):
        num[v] = order
        order += 1
        low[v] = num[v]
        n_child = 0
        for u in VtoV[v]:
            if u == par:
                continue
            if num[u] != -1:
                low[v] = min(low[v], num[u])
                continue
            n_child += 1
            order = dfs(u, v, order)
            low[v] = min(low[v], low[u])
            if par != -1 and num[v] <= low[u]:
                cutvertices.append(v)
        if par == -1 and n_child >= 2:
            cutvertices.append(v)
        return order
    dfs(0, -1, 0)
    return cutvertices


def deleteUnnecessaryVertices(N: int, edges: list, degrees: list, show: list, multiedges: list, is_remain: bool):
    diagonal = calcDiagonalVertices(N, edges)
    # show = True の頂点で構成されるグラフを作る
    VtoV_deleted = [[] for _ in range(N)]
    for e in edges:
        if not show[e[0]] or not show[e[1]]:
            continue
        if any([show[v] for v in diagonal[e]]):
            VtoV_deleted[e[0]].append(e[1])
    cutvertices = getCutVertices(N, VtoV_deleted)
    # カット点に隣接している点とそれに接続する辺は show が False でも残す。
    cutvertices_neighbors = []
    for v in range(N):
        for c in cutvertices:
            if not show[v] and (v, c) in edges:
                cutvertices_neighbors.append(v)
    for e in edges:
        if (e[0] in cutvertices_neighbors and (e[1] in cutvertices or any([u in cutvertices for u in diagonal[e]]))) or \
           (e[1] in cutvertices_neighbors and (e[0] in cutvertices or any([u in cutvertices for u in diagonal[e]]))):
            VtoV_deleted[e[0]].append(e[1])
    
    # いらない点を削除したグラフを作成する。
    new_N = 0

    new_index = [-1 for _ in range(N)]
    for v in range(N):
        if show[v] or v in cutvertices_neighbors:
            new_index[v] = new_N
            new_N += 1
    # is_remain == True (reducibility のチェックにいらない頂点も残す) のときは、
    # cartwheel の中心を区別するために 0 番目の頂点(中心) の番号は変わらないでほしい。
    assert not is_remain or new_index[0] == 0 

    new_edges = []
    for v in range(N):
        for u in VtoV_deleted[v]:
            if new_index[u] >= 0 and new_index[v] >= 0:
                new_edges.append((new_index[u], new_index[v]))

    new_degrees = [0 for _ in range(new_N)]
    for v in range(N):
        if show[v]:
            new_degrees[new_index[v]] = degrees[v]
        elif v in cutvertices_neighbors:
            new_degrees[new_index[v]] = INF
    
    new_multiedges = []
    for e in multiedges:
        if not show[e[0]] or not show[e[1]]:
            continue
        if any([show[v] for v in diagonal[e]]):
            e0, e1 = new_index[e[0]], new_index[e[1]]
            assert e0 >= 0 and e1 >= 0
            new_multiedges.append((min(e0, e1), max(e0, e1)))
    
    return new_N, new_edges, new_degrees, new_multiedges
    

# 
def getCycleFromNeighbors(N: int, cycle_neighbors: list):
    si = -1
    count = 0
    for i in range(N):
        assert len(cycle_neighbors[i]) == 0 or len(cycle_neighbors[i]) == 2, f"{len(cycle_neighbors[i])}"
        if len(cycle_neighbors[i]) == 2:
            si = i
            count += 1
    assert si != -1
    outer_face = []
    outer_face.append(si)
    i = cycle_neighbors[si][0]
    while i != si:
        new_i = -1
        if cycle_neighbors[i][0] == outer_face[-1]:
            new_i = cycle_neighbors[i][1]
        elif cycle_neighbors[i][1] == outer_face[-1]:
            new_i = cycle_neighbors[i][0]
        assert new_i != -1
        outer_face.append(i)
        i = new_i
    assert outer_face[-1] == cycle_neighbors[si][1]
    assert len(outer_face) == count
    return outer_face


# outer_face を得る。 
def calcOuterFace(N: int, edges: list):
    diagonal = calcDiagonalVertices(N, edges)
    cycle_neighbors = [[] for _ in range(N)]
    for e in edges:
        if len(diagonal[e]) == 1 and e[0] < e[1]:
            cycle_neighbors[e[0]].append(e[1])
            cycle_neighbors[e[1]].append(e[0])
    return getCycleFromNeighbors(N, cycle_neighbors)


# outer_face の頂点が残り何本の辺を outer_face に出すかを計算する
def calcNumNeighbors(N: int, edges: list, degrees: list, outer_face: list, multiedges: list):
    VtoV_num = [0 for _ in range(N)]
    incident_multi = [0 for _ in range(N)]
    for e in edges:
        if e[0] < e[1]:
            VtoV_num[e[0]] += 1
            VtoV_num[e[1]] += 1
    for e in multiedges:
        assert incident_multi[e[0]] == 0 and incident_multi[e[1]] == 0
        incident_multi[e[0]] = 1
        incident_multi[e[1]] = 1
    
    num_neighbors = [0 for _ in range(len(outer_face))]
    for i in range(len(outer_face)):
        if degrees[outer_face[i]] == INF: # カット点に隣接している ring の頂点
            num_neighbors[i] = INF
        else:
            num_neighbors[i] = degrees[outer_face[i]] - VtoV_num[outer_face[i]] - incident_multi[outer_face[i]]
            assert num_neighbors[i] >= 0, f"num_neighbors[i]: {num_neighbors[i]}, degrees[outer_face[i]]: {degrees[outer_face[i]]}, VtoV_num[outer_face[i]]: {VtoV_num[outer_face[i]]} incident_multi[outer_face[i]]: {incident_multi[outer_face[i]]}"
            if i > 0:
                assert num_neighbors[i] > 0 or num_neighbors[i - 1] > 0
            if i == len(outer_face) - 1:
                assert num_neighbors[i] > 0 or num_neighbors[0] > 0
    return num_neighbors


def handleZeroNeighbor(i: int, l: int, outer_face: list, num_new_neighbors: list, 
                       is_outer_face: list, before: int, after: int, edges: list):
    assert num_new_neighbors[i] == 0 and is_outer_face[i]
    is_outer_face[i] = False
    edges.append((outer_face[before[i]], outer_face[after[i]]))
    edges.append((outer_face[after[i]], outer_face[before[i]]))
    num_new_neighbors[before[i]] -= 1
    num_new_neighbors[after[i]] -= 1
    before[after[i]] = before[i]
    after[before[i]] = after[i]
    assert is_outer_face[before[i]] and is_outer_face[after[i]]
    if num_new_neighbors[before[i]] == 0:
        handleZeroNeighbor(before[i], l, outer_face, num_new_neighbors, is_outer_face, before, after, edges)
    if is_outer_face[after[i]] and num_new_neighbors[after[i]] == 0:
        handleZeroNeighbor(after[i], l, outer_face, num_new_neighbors, is_outer_face, before, after, edges)
    return


# ring を補完する
def complementRing(N: int, edges: list, degrees: list, multiedges: list):
    outer_face = calcOuterFace(N, edges)
    num_new_neighbors = calcNumNeighbors(N, edges, degrees, outer_face, multiedges)

    # num_new_neighbors[v] = 0 の頂点の処理
    l = len(outer_face)
    is_outer_face = [True for _ in range(l)]
    before = [(i + l - 1) % l for i in range(l)]
    after = [(i + 1) % l for i in range(l)]
    for i in range(l):
        if is_outer_face[i] and num_new_neighbors[i] == 0:
            handleZeroNeighbor(i, l, outer_face, num_new_neighbors, is_outer_face, before, after, edges)
    outer_face2 = []
    num_new_neighbors2 = []
    for i in range(l):
        if is_outer_face[i]:
            outer_face2.append(outer_face[i])
            num_new_neighbors2.append(num_new_neighbors[i])
    outer_face, num_new_neighbors = outer_face2, num_new_neighbors2

    l = len(outer_face)
    n_group = 0
    left_group = [-1 for _ in range(l)]
    right_group = [-1 for _ in range(l)]

    # assert num_new_neighbors[0] < INF - 5 or num_new_neighbors[l-1] < INF - 5
    # が満たされるように outer_face を rotate する。
    offset = 0
    while offset < l:
        if num_new_neighbors[offset] < INF - 5 and num_new_neighbors[(offset + l - 1) % l] < INF - 5:
            break
        offset += 1
    assert offset < l
    outer_face = outer_face[offset:] + outer_face[:offset]
    num_new_neighbors = num_new_neighbors[offset:] + num_new_neighbors[:offset]
    assert num_new_neighbors[0] < INF - 5 or num_new_neighbors[l-1] < INF - 5

    i = 0
    while i < l:
        assert num_new_neighbors[i] > 0
        if num_new_neighbors[i] == 1:
            left_group[i] = n_group
            right_group[i] = n_group
        elif num_new_neighbors[i] >= INF - 5:
            left_group[i] = n_group
            while num_new_neighbors[(i + 1) % l] >= INF - 5:
                i += 1
            n_group += 1
            right_group[i] = n_group
        else:
            left_group[i] = n_group
            n_group += 1
            right_group[i] = n_group
        i += 1
    assert n_group > 0
    i = l - 1
    while num_new_neighbors[i] == 1:
        right_group[i] = left_group[0]
        left_group[i] = left_group[0]
        i -= 1
    right_group[i] = left_group[0]

    group_neighbors = [i for i in range(N, N + n_group)]

    new_N = N + n_group
    new_edges = edges.copy()
    new_degrees = degrees.copy() + [INF for _ in range(n_group)]
    new_multiedges = multiedges.copy()
    def new_vertex():
        nonlocal new_N
        v = new_N
        new_N += 1
        new_degrees.append(INF)
        return v
    def add_edge(u: int, v: int):
        new_edges.append((u, v))
        new_edges.append((v, u))
        return
    for i in range(l):
        if num_new_neighbors[i] == 1:
            add_edge(outer_face[i], group_neighbors[left_group[i]])
            assert left_group[i] == right_group[i]
        else:
            if left_group[i] != -1:
                first = group_neighbors[left_group[i]]
                add_edge(outer_face[i], first)
            if right_group[i] != -1:
                last = group_neighbors[right_group[i]]
                add_edge(outer_face[i], last)
            if num_new_neighbors[i] >= INF - 5:
                continue
            for _ in range(num_new_neighbors[i] - 2):
                v = new_vertex()
                add_edge(v, first)
                add_edge(v, outer_face[i])
                first = v
            add_edge(first, last)
    
    assert new_N == len(new_degrees)
    return new_N, new_edges, new_degrees, new_multiedges


# conf のインデックスを付け直す。
def genConf(N: int, edges: list, degrees: list, multiedges: list, is_remain: bool):
    # ring
    ring_neighbors = [[] for _ in range(N)]
    for e in edges:
        if degrees[e[0]] >= INF - 5 and degrees[e[1]] >= INF - 5:
            ring_neighbors[e[0]].append(e[1])
    ring = getCycleFromNeighbors(N, ring_neighbors)
    R = len(ring)

    new_index = [-1 for _ in range(N)]
    idx = 0
    for v in ring:
        new_index[v] = idx
        idx += 1
    for v in range(N):
        if new_index[v] == -1:
            new_index[v] = idx
            idx += 1
    assert idx == N
    # is_remain == True (reducibility のチェックにいらない頂点も残す) のときは、
    # cartwheel の中心を区別するために 0 番目の頂点(中心) の番号は R になって欲しい
    assert not is_remain or new_index[0] == R

    new_VtoV = [[] for _ in range(N)]
    for e in edges:
        e0, e1 = new_index[e[0]], new_index[e[1]]
        new_VtoV[e0].append(e1)

    new_multiedges = []
    incident_multi = [0 for _ in range(N)]
    for e in multiedges:
        e0, e1 = new_index[e[0]], new_index[e[1]]
        assert e0 >= R and e1 >= R
        new_multiedges.append((min(e0, e1), max(e0, e1)))
        incident_multi[e0] += 1
        incident_multi[e1] += 1
        assert incident_multi[e0] == 1 and incident_multi[e1] == 1

    for v in range(N):
        if new_index[v] >= R:
            assert degrees[v] == len(new_VtoV[new_index[v]]) + incident_multi[new_index[v]]
    
    return N, R, new_VtoV, new_multiedges


# テキストにする
def genConTxt(N: int, R: int, VtoV: list, multiedges: list):
    text = f"\n{N} {R}\n"
    for v in range(R, N):
        text_v = f"{v+1} {len(VtoV[v])}"
        for u in VtoV[v]:
            text_v += f" {u + 1}"
        text += text_v + "\n"
    if len(multiedges) > 0:
        text += f"{len(multiedges)}\n"
        for e in multiedges:
            text += f"{e[0]+1} {e[1]+1}\n"
    return text


def cartWheel2Conf(cartwheel_str: str, is_remain: bool, is_parallel: bool):
    # cartwheel の生成
    N, nodes, edges, degrees, show, multiedges = cartWheel(cartwheel_str, is_parallel)
    # 外平面に 4 本以上の辺を出している頂点 v は show[v] = False とする。
    markUnnecessaryVertices(N, edges, degrees, show, multiedges, is_remain)
    # いらない頂点 (show[v] = False) を消す。
    # ただし、カット点がある場合はカット点に隣接している頂点のみ残す。
    # このときそのような頂点と show[v] = True となる頂点が多重辺で結ばれている場合は
    # その多重辺を通常の辺に直す。
    N, edges, degrees, multiedges = deleteUnnecessaryVertices(N, edges, degrees, show, multiedges, is_remain)
    # ring を補完する 
    # degrees[v] が INF になっている頂点が ring の頂点
    N, edges, degrees, multiedges = complementRing(N, edges, degrees, multiedges)
    # index をつけ直す。
    N, R, VtoV, multiedges = genConf(N, edges, degrees, multiedges, is_remain)
    text = genConTxt(N, R, VtoV, multiedges)
    return text


def confFromLog(file: str, is_remain: bool, is_parallel: bool):
    pattern = r'cartwheel : '
    cPattern = r'charge \(initial, receive, send, result\) : -?(\d+), (\d+), (\d+), -?(\d+)'
    prog = re.compile(pattern)
    cProg = re.compile(cPattern)
    conf = None
    cartwheels = []
    with open(file) as f:
        for line in f:
            cResult = cProg.search(line)
            if cResult is not None:
                charge = int(cResult[4])
                if charge == 0:
                    cartwheels.append(conf)
            result = prog.search(line)
            if result is None:
                continue
            conf = cartWheel2Conf(line[result.end(0):], is_remain, is_parallel) 
    return cartwheels 


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('inputPath', help="A file that contains log")
    parser.add_argument('outputDir')
    parser.add_argument('--prefix', default="conf")
    parser.add_argument('--index', type=int, default=0)
    parser.add_argument('--remain', action='store_true', help="A flag to decide whether unnecessary vertices for reducible check remain")
    parser.add_argument('-p', '--parallel', action='store_true')
    args = parser.parse_args()
    if not os.path.isdir(args.outputDir):
        print("output Dir must be directory!")
        exit(1)
    confs = confFromLog(args.inputPath, args.remain, args.parallel)
    for i, conf in enumerate(confs, int(args.index)):
        path = os.path.join(args.outputDir, f"{args.prefix}{i:04}.conf")
        with open(path, 'w') as f:
            f.write(conf)
        
