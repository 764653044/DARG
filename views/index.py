from django.shortcuts import render, redirect
from django import forms
from drug_addiction.models import UploadedFile
from django.http import JsonResponse, HttpResponse
import networkx as nx
import io
import numpy as np
import igraph as ig
import csv
import pymysql
from operator import itemgetter

config = {
    'host': '127.0.0.1',
    'port': 3306,
    'user': 'root',
    'password': '123456',
    'database': 'database',
    'charset': 'utf8',
}

def index(request):
    return render(request, 'index.html')

def read_uploaded_file(file):
    # 使用TextIOWrapper将文件流转换为文本文件对象，指定编码为utf-8
    text_file = io.TextIOWrapper(file, encoding='utf-8')
    # 逐行读取文件内容
    lines = []
    for line in text_file:
        lines.append(line.strip())  # 去除每行末尾的换行符
    return lines

def upload(request:"HttpRequest"):
    if request.method == 'POST':
        species = request.POST.get('species')
        genename = request.POST.get('gene')
        threshold = request.POST.get('threshold')
        request.session['threshold'] = threshold
        genename_list = genename.split()
        request.session['user_update_data'] = genename_list
        request.session['species'] = species
        file = request.FILES.get("fileUpLoad")  # 获取文件
        if file:
            user_update_data = read_uploaded_file(file)
            # 现在，'lines' 是包含文件每一行内容的列表，您可以对其进行进一步处理
            request.session['user_update_data'] = user_update_data
        return redirect('success')
    return render(request, 'upload.html')

def radiality(graph, mode="all"):
    sp = dict(nx.all_pairs_shortest_path_length(graph))
    n = len(graph)
    diam = max(max(val.values()) for val in sp.values())
    radiality_dict = {}
    for node in graph.nodes():
        neighbors = list(graph.neighbors(node))
        if len(neighbors) == 0:
            radiality_dict[node] = 0
        else:
            x = []
            for neighbor in neighbors:
                if neighbor in sp[node]:
                    x.extend([diam + 1 - sp[node][neighbor]])
            radiality_dict[node] = sum(x) / (n - 1)
    return radiality_dict

def epc_centrality(graph, threshold=0.5, iterations=1000, random_seed=None):
    if random_seed is not None:
        np.random.seed(random_seed)
    epc_scores = np.zeros(graph.vcount())  # 初始化 EPC 中心性得分为0
    for i in range(iterations):
        g = graph.copy()
        g.es['weight'] = np.random.uniform(0, 1, size=g.ecount())
        edges_to_remove = [edge.index for edge in g.es if edge['weight'] < threshold]
        g.delete_edges(edges_to_remove)
        components = g.clusters()
        membership = {v: c for c, component in enumerate(components) for v in component}
        epc = np.zeros(graph.vcount())
        for v in g.vs:
            for t in g.vs:
                if membership[v.index] == membership[t.index]:
                    epc[v.index] += 1
        epc_scores += epc
    epc_scores /= iterations
    return dict(zip(graph.vs['name'], epc_scores))

def semilocal(graph, vids=None, mode="all"):
    if vids is None:
        vids = graph.nodes()  # 使用所有节点作为默认值
    res = {}
    for v in vids:
        vNeighbors = list(graph.neighbors(v))
        sl = 0
        for vv in vNeighbors:
            vvNeighbors = list(graph.neighbors(vv))
            for vvv in vvNeighbors:
                sl += len(list(graph.neighbors(vvv)))
        res[v] = sl
    return res

def laplacian(graph, mode="all", loops=True):
    res = []
    for aNode in graph.vs:
        aNeighbors = aNode.neighbors(mode=mode)
        deg = sum(aNeighbor.degree(mode=mode, loops=loops) for aNeighbor in aNeighbors)
        aDeg = aNode.degree(mode=mode, loops=loops)
        res.append(aDeg**2 + aDeg + 2 * deg)
    node_names = graph.vs["name"]
    result_dict = dict(zip(node_names, res))
    return result_dict

def success_view(request):
    user_update_data = request.session.get('user_update_data')
    species = request.session.get('species')
    threshold = request.session.get('threshold')
    user_update_data_tuple = tuple(user_update_data)

    db = pymysql.connect(**config)
    cursor = db.cursor()
    if species == 'Human':
        sql = "SELECT gene1, gene2 FROM `ppi-human-links-v12.0` WHERE gene1 IN " + str(user_update_data_tuple) + " AND gene2 IN " + str(
            user_update_data_tuple) + "AND combined_score >=" + str(threshold)
    else:
        sql = "SELECT gene1, gene2 FROM `ppi-mouse-links-v12.0` WHERE gene1 IN " + str(
            user_update_data_tuple) + " AND gene2 IN " + str(
            user_update_data_tuple) + "AND combined_score >=" + str(threshold)
    cursor.execute(sql)
    data = cursor.fetchall()
    db.close()

    ppi_data1 = []
    ppi_data2 = []
    for j in data:
        ppi_data1.append(j[0])
        ppi_data2.append(j[1])
    zipped_ppi = zip(ppi_data1, ppi_data2)

    user_update_data0 = data
    # user_update_data0 = [item.split(',') for item in user_update_data]

    G = nx.Graph()
    # for edge in user_update_data0:
    #     G.add_edge(edge[0], edge[1])
    G.add_edges_from(user_update_data0)
# Degree
    node_degrees = {}
    for edge in user_update_data0:
        source_node = edge[0]
        target_node = edge[1]
        if source_node not in node_degrees:
            node_degrees[source_node] = 0
        if target_node not in node_degrees:
            node_degrees[target_node] = 0
        node_degrees[source_node] += 1
        node_degrees[target_node] += 1
# MNC
    mnc_centralities = {}
    for node in G.nodes():
        mnc = len(max(nx.connected_components(G.subgraph(G.neighbors(node))), key=len))
        mnc_centralities[node] = mnc
# Katz
    katz_centrality = nx.katz_centrality(G, alpha=0.01)  # alpha是衡量路径长度的参数，可以根据需要调整
# Radiality
    radiality_dict = radiality(G, mode="in")
# Semilocal
    semilocal_scores = semilocal(G, mode="in")
# Laplacian
    G = ig.Graph.TupleList(user_update_data0, directed=False)
    laplacian_centrality = laplacian(G, mode="in")
# EPC
    epc_results = epc_centrality(G, threshold=0.5, random_seed=42)

    centrality_number = len(semilocal_scores) // 10

    semilocal_dict = dict(sorted(semilocal_scores.items()))
    semilocal_dict_values = semilocal_dict.values()
    gene_names = semilocal_dict.keys()

    degree_dict = dict(sorted(node_degrees.items()))
    degree_dict_values = degree_dict.values()

    mnc_dict = dict(sorted(mnc_centralities.items()))
    mnc_dict_values = mnc_dict.values()

    katz_dict = dict(sorted(katz_centrality.items()))
    katz_dict_values = katz_dict.values()

    radiality_dict = dict(sorted(radiality_dict.items()))
    radiality_dict_values = radiality_dict.values()

    laplacian_dict = dict(sorted(laplacian_centrality.items()))
    laplacian_dict_values = laplacian_dict.values()

    epc_dict = dict(sorted(epc_results.items()))
    epc_dict_values = epc_dict.values()

    result = []
    # 使用zip函数将两个dict_values合并
    for val1, val2, val3, val4, val5, val6, val7, val8 in zip(gene_names, degree_dict_values, radiality_dict_values, epc_dict_values, mnc_dict_values,
                                katz_dict_values, semilocal_dict_values, laplacian_dict_values):
        result.append([val1, val2, val3, val4, val5, val6, val7, val8])

    semilocal_list = list(semilocal_scores.items())
    degree_list = list(node_degrees.items())
    mnc_list = list(mnc_centralities.items())
    katz_list = list(katz_centrality.items())
    radiality_list = list(radiality_dict.items())
    laplacian_list = list(laplacian_centrality.items())
    epc_list = list(epc_results.items())

    # 对每个列表按照元组的第一个元素进行排序
    semilocal_list.sort(key=lambda x: x[1], reverse=True)
    semilocal_list = semilocal_list[:centrality_number]
    semilocal_top_node = []
    for i in semilocal_list:
        semilocal_top_node.append(i[0])
    set_semilocal = set(semilocal_top_node)
    semilocal_list = [[item[0], item[1]] for item in semilocal_list]

    degree_list.sort(key=lambda x: x[1], reverse=True)
    degree_list = degree_list[:centrality_number]
    degree_top_node = []
    for i in degree_list:
        degree_top_node.append(i[0])
    set_degree = set(degree_top_node)
    degree_list = [[item[0], item[1]] for item in degree_list]

    mnc_list.sort(key=lambda x: x[1], reverse=True)
    mnc_list = mnc_list[:centrality_number]
    mnc_top_node = []
    for i in mnc_list:
        mnc_top_node.append(i[0])
    set_mnc = set(mnc_top_node)
    mnc_list = [[item[0], item[1]] for item in mnc_list]

    katz_list.sort(key=lambda x: x[1], reverse=True)
    katz_list = katz_list[:centrality_number]
    katz_top_node = []
    for i in katz_list:
        katz_top_node.append(i[0])
    set_katz = set(katz_top_node)
    katz_list = [[item[0], item[1]] for item in katz_list]

    radiality_list.sort(key=lambda x: x[1], reverse=True)
    radiality_list = radiality_list[:centrality_number]
    radiality_top_node = []
    for i in radiality_list:
        radiality_top_node.append(i[0])
    set_radiality = set(radiality_top_node)
    radiality_list = [[item[0], item[1]] for item in radiality_list]

    laplacian_list.sort(key=lambda x: x[1], reverse=True)
    laplacian_list = laplacian_list[:centrality_number]
    laplacian_top_node = []
    for i in laplacian_list:
        laplacian_top_node.append(i[0])
    set_laplacian = set(laplacian_top_node)
    laplacian_list = [[item[0], item[1]] for item in laplacian_list]

    epc_list.sort(key=lambda x: x[1], reverse=True)
    epc_list = epc_list[:centrality_number]
    epc_top_node = []
    for i in epc_list:
        epc_top_node.append(i[0])
    set_epc = set(epc_top_node)
    epc_list = [[item[0], item[1]] for item in epc_list]

    intersection = set_degree & set_radiality & set_epc & set_mnc & set_katz & set_laplacian & set_semilocal
    intersection_list = []
    for i in intersection:
        intersection_list.append(i)

    for i in result:
        if i[0] in intersection_list:
            i.append('Key gene')
    key_gene_items = [item for item in result if item[-1] == 'Key gene']
    not_key_list = [item for item in result if item[-1] != 'Key gene']
    all_centrality_list = key_gene_items + not_key_list
    request.session['all_centrality_list'] = all_centrality_list

    key_degree_list = [item for item in degree_list if item[0] in intersection_list]
    key_radiality_list = [item for item in radiality_list if item[0] in intersection_list]
    key_epc_list = [item for item in epc_list if item[0] in intersection_list]
    key_mnc_list = [item for item in mnc_list if item[0] in intersection_list]
    key_katz_list = [item for item in katz_list if item[0] in intersection_list]
    key_semilocal_list = [item for item in semilocal_list if item[0] in intersection_list]
    key_laplacian_list = [item for item in laplacian_list if item[0] in intersection_list]

    key_gene = []
    key_degree = []
    for i in key_degree_list:
        key_gene.append(i[0])
        key_degree.append(i[1])
    key_radiality = []
    for i in key_radiality_list:
        key_radiality.append(i[1])
    key_epc = []
    for i in key_epc_list:
        key_epc.append(i[1])
    key_mnc = []
    for i in key_mnc_list:
        key_mnc.append(i[1])
    key_katz = []
    for i in key_katz_list:
        key_katz.append(i[1])
    key_semilocal = []
    for i in key_semilocal_list:
        key_semilocal.append(i[1])
    key_laplacian = []
    for i in key_laplacian_list:
        key_laplacian.append(i[1])
    zipped_data = zip(key_gene, key_degree, key_radiality, key_epc, key_mnc, key_katz, key_semilocal, key_laplacian)
    return render(request, 'analyseResult.html', locals())

def download(request):
    all_centrality_list = request.session.get('all_centrality_list')
    filename = 'all_gene_centrality.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)
    row = []
    writer.writerow(
        ['Gene', 'Degree', 'Radiality', 'EPC', 'MNC', 'Katz', 'Semilocal', 'Laplacian', 'Key'])
    for i in all_centrality_list:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def keygene(request):
    gene = request.GET.get('gene')
    # 获取同源基因
    if any(char.islower() for char in gene):
        db = pymysql.connect(**config)
        cursor4 = db.cursor()
        sql4 = "SELECT human,mouse_ID,human_ID FROM homologene WHERE mouse = " + "'" + str(gene) + "'"
        cursor4.execute(sql4)
        data4 = cursor4.fetchall()
        if data4 != ():
            homologene_human = data4[0][0]
            mouse_ID = data4[0][1]
            human_ID = data4[0][2]
        else:
            homologene_human = gene
        db.close()
    else:
        db = pymysql.connect(**config)
        cursor8 = db.cursor()
        sql8 = "SELECT mouse,mouse_ID,human_ID FROM homologene WHERE human = " + "'" + str(gene) + "'"
        cursor8.execute(sql8)
        data8 = cursor8.fetchall()
        if data8 != ():
            homologene_mouse = data8[0][0]
            mouse_ID = data8[0][1]
            human_ID = data8[0][2]
        else:
            homologene_mouse = gene
        homologene_human = gene
        db.close()

    # 获取ctd-cocaine
    cocaine_list = []
    db = pymysql.connect(**config)
    cursor5 = db.cursor()
    sql5 = "SELECT `Inference Score`,`Reference Count` FROM ctd_cocaine WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor5.execute(sql5)
    data5 = cursor5.fetchall()
    for i in data5:
        cocaine_list.append(float(i[0]))
        cocaine_list.append(float(i[1]))
    db.close()

    # 获取ctd-heroin
    heroin_list = []
    db = pymysql.connect(**config)
    cursor6 = db.cursor()
    sql6 = "SELECT `Inference Score`,`Reference Count` FROM ctd_heroin WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor6.execute(sql6)
    data6 = cursor6.fetchall()
    for i in data6:
        heroin_list.append(float(i[0]))
        heroin_list.append(float(i[1]))
    db.close()

    # 获取ctd-morphine
    morphine_list = []
    db = pymysql.connect(**config)
    cursor7 = db.cursor()
    sql7 = "SELECT `Inference Score`,`Reference Count` FROM ctd_morphine WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor7.execute(sql7)
    data7 = cursor7.fetchall()
    for i in data7:
        morphine_list.append(float(i[0]))
        morphine_list.append(float(i[1]))
    db.close()

    # 获取ctd-amphetamine
    amphetamine_list = []
    db = pymysql.connect(**config)
    cursor14 = db.cursor()
    sql14 = "SELECT `Inference Score`,`Reference Count` FROM ctd_amphetamine WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor14.execute(sql14)
    data14 = cursor14.fetchall()
    for i in data14:
        amphetamine_list.append(float(i[0]))
        amphetamine_list.append(float(i[1]))
    db.close()

    # 获取ctd-cannabis
    cannabis_list = []
    db = pymysql.connect(**config)
    cursor13 = db.cursor()
    sql13 = "SELECT `Inference Score`,`Reference Count` FROM ctd_cannabis WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor13.execute(sql13)
    data13 = cursor13.fetchall()
    for i in data13:
        cannabis_list.append(float(i[0]))
        cannabis_list.append(float(i[1]))
    db.close()

    # 获取ctd-ethanol
    ethanol_list = []
    db = pymysql.connect(**config)
    cursor12 = db.cursor()
    sql12 = "SELECT `Inference Score`,`Reference Count` FROM ctd_ethanol WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor12.execute(sql12)
    data12 = cursor12.fetchall()
    for i in data12:
        ethanol_list.append(float(i[0]))
        ethanol_list.append(float(i[1]))
    db.close()

    # 获取ctd-nicotine
    nicotine_list = []
    db = pymysql.connect(**config)
    cursor11 = db.cursor()
    sql11 = "SELECT `Inference Score`,`Reference Count` FROM ctd_nicotine WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor11.execute(sql11)
    data11 = cursor11.fetchall()
    for i in data11:
        nicotine_list.append(float(i[0]))
        nicotine_list.append(float(i[1]))
    db.close()

    # 获取ctd中药物与基因的文献信息
    db = pymysql.connect(**config)
    cursor8 = db.cursor()
    sql8 = "SELECT ChemicalName,ChemicalID,Organism,Interaction,InteractionActions,PubMedIDs FROM ctd_drug_gene_link WHERE `GeneSymbol` = " + "'" + str(
        homologene_human) + "'"
    cursor8.execute(sql8)
    data8 = cursor8.fetchall()
    db.close()

    """ 获取hpa数据库中基因的脑区表达量 """
    # 人-HPA
    db = pymysql.connect(**config)
    cursor15 = db.cursor()
    sql15 = "SELECT `data` FROM hpa_rna_human_brain WHERE `gene` =" + "'" + str(gene) + "'"
    cursor15.execute(sql15)
    data15 = cursor15.fetchall()
    db.close()
    if data15:
        human_hpa_list = data15[0][0]

    # 人-HPA-热图
    db = pymysql.connect(**config)
    cursor20 = db.cursor()
    sql20 = "SELECT `data` FROM hpa_rna_human_brain_heatmap WHERE `gene` =" + "'" + str(gene) + "'"
    cursor20.execute(sql20)
    data20 = cursor20.fetchall()
    db.close()
    if data20:
        human_hpa_heatmap_list = data20[0][0]

    # 人-GTEX
    db = pymysql.connect(**config)
    cursor18 = db.cursor()
    sql18 = "SELECT `data` FROM hpa_rna_human_brain_gtex WHERE `gene` =" + "'" + str(gene) + "'"
    cursor18.execute(sql18)
    data18 = cursor18.fetchall()
    db.close()
    if data18:
        human_gtex_list = data18[0][0]

    # 人-GTEX-热图
    db = pymysql.connect(**config)
    cursor21 = db.cursor()
    sql21 = "SELECT `data` FROM hpa_rna_human_brain_gtex_heatmap WHERE `gene` =" + "'" + str(gene) + "'"
    cursor21.execute(sql21)
    data21 = cursor21.fetchall()
    db.close()
    if data21:
        human_gtex_heatmap_list = data21[0][0]

    # 人-FANTOM
    db = pymysql.connect(**config)
    cursor19 = db.cursor()
    sql19 = "SELECT `data` FROM hpa_rna_human_brain_fantom WHERE `gene` =" + "'" + str(gene) + "'"
    cursor19.execute(sql19)
    data19 = cursor19.fetchall()
    db.close()
    if data19:
        human_fantom_list = data19[0][0]

    # 人-FANTOM-热图
    db = pymysql.connect(**config)
    cursor22 = db.cursor()
    sql22 = "SELECT `data` FROM hpa_rna_human_brain_fantom_heatmap WHERE `gene` =" + "'" + str(gene) + "'"
    cursor22.execute(sql22)
    data22 = cursor22.fetchall()
    db.close()
    if data22:
        human_fantom_heatmap_list = data22[0][0]

    # 鼠
    db = pymysql.connect(**config)
    cursor16 = db.cursor()
    sql16 = "SELECT `data` FROM hpa_rna_mouse_brain WHERE `gene` =" + "'" + str(gene) + "'"
    cursor16.execute(sql16)
    data16 = cursor16.fetchall()
    db.close()
    if data16:
        mouse_hpa_list = data16[0][0]

    # 鼠-热图
    db = pymysql.connect(**config)
    cursor23 = db.cursor()
    sql23 = "SELECT * FROM hpa_rna_mouse_brain_heatmap WHERE `gene` =" + "'" + str(gene) + "'"
    cursor23.execute(sql23)
    data23 = cursor23.fetchall()
    db.close()
    if data23:
        mouse_hpa_heatmap_list = data23[0]
        # 第一个大图
        mouse_hpa_heatmap_list1 = mouse_hpa_heatmap_list[1:]
        mouse_hpa_heatmap_list1 = str(mouse_hpa_heatmap_list1).replace("'", "")
        mouse_hpa_heatmap_list1 = str(mouse_hpa_heatmap_list1)[1:-1]

    # 鼠-allen
    db = pymysql.connect(**config)
    cursor17 = db.cursor()
    sql17 = "SELECT `data` FROM hpa_rna_mouse_brain_allen WHERE `gene` =" + "'" + str(gene) + "'"
    cursor17.execute(sql17)
    data17 = cursor17.fetchall()
    db.close()
    if data17:
        mouse_allen_list = data17[0][0]

    # 鼠-allen-热图
    db = pymysql.connect(**config)
    cursor24 = db.cursor()
    sql24 = "SELECT * FROM hpa_rna_mouse_brain_allen_heatmap WHERE `gene` =" + "'" + str(gene) + "'"
    cursor24.execute(sql24)
    data24 = cursor24.fetchall()
    db.close()
    if data24:
        mouse_allen_heatmap_list = data24[0]
        # 第一个大图
        mouse_allen_heatmap_list1 = mouse_allen_heatmap_list[1:]
        mouse_allen_heatmap_list1 = str(mouse_allen_heatmap_list1).replace("'", "")
        mouse_allen_heatmap_list1 = str(mouse_allen_heatmap_list1)[1:-1]

    return render(request, 'key_gene.html', locals())