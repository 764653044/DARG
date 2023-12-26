from django.shortcuts import render, redirect
from django.http import JsonResponse, HttpResponse
from drug_addiction.models import SearchCriteria
from django.core import serializers
import csv
import pymysql

config = {
    'host': '127.0.0.1',
    'port': 3306,
    'user': 'root',
    'password': '123456',
    'database': 'database',
    'charset': 'utf8',
}

def getInfo(request):
    pid = request.GET.get('pid', -1)
    pid = int(pid)
    areaList = SearchCriteria.objects.filter(parentid=pid)
    jareaList = serializers.serialize('json', areaList)
    return JsonResponse({'jareaList': jareaList})

def help(request):

    return render(request, 'help.html')

def download(request):
    a = request.GET.get('a')
    b = request.GET.get('b')
    c = request.GET.get('c')
    d = request.GET.get('d')
    p = request.GET.get('p')
    fc = request.GET.get('fc')
    if not p:
        p = '0.05'
        fc = '1.2'
    request.session['p'] = p
    request.session['fc'] = fc



    table_name = str(a) + "_" + str(b) + "_" + str(c) + "_" + str(d)
    table_name_lower = table_name.lower()
    if a == None:
        table_name = "Human_Cocaine_Midbrain_GSE54839"
    request.session['table_name'] = table_name

    db = pymysql.connect(**config)
    cursor = db.cursor()
    sql = "SHOW TABLES LIKE '"+str(a)+"_"+str(b)+"_"+str(c)+"_"+str(d)+"_ppi%'"
    cursor.execute(sql)
    data = cursor.fetchall()
    db.close()
    data_list = []
    for i in data:
        data_list.append(i[0])

    return render(request, 'download.html', locals())

def methylationdownload(request):
    a = request.GET.get('a')
    b = request.GET.get('b')
    c = request.GET.get('c')
    p = request.GET.get('p')
    deltaBeta = request.GET.get('deltaBeta')
    if not p:
        p = '0.05'
        deltaBeta = '0.05'
    request.session['p'] = p
    request.session['deltaBeta'] = deltaBeta

    table_name = str(a) + "_" + str(b) + "_" + str(c)
    table_name_lower = table_name.lower()
    if a == None:
        table_name = "Cocaine_PeripheralBlood_GSE77056"
    request.session['table_name'] = table_name

    db = pymysql.connect(**config)
    cursor = db.cursor()
    sql = "SHOW TABLES LIKE 'm_" + str(a) + "_" + str(b) + "_" + str(c) + "_go%'"
    cursor.execute(sql)
    data = cursor.fetchall()
    db.close()
    data_list = []
    for i in data:
        data_list.append(i[0])

    return render(request, 'download_methylation.html', locals())

def download1(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    fc = request.session.get('fc')

    db = pymysql.connect(**config)
    cursor0 = db.cursor()
    sql0 = "SELECT type FROM dataset_type WHERE table_name = " + "'" + str(table_name) + "'"
    cursor0.execute(sql0)
    type_dataset = cursor0.fetchall()
    type_dataset = type_dataset[0][0]
    db.close()

    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if type_dataset == 'limma':
        sql1 = "SELECT * FROM " + str(table_name) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ")"
    else:
        sql1 = "SELECT * FROM " + str(table_name) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ")"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    if type_dataset == 'limma':
        writer.writerow(['Gene', 'logFC', 'AveExpr', 't', 'P.Value',
                         'adj.P.Val', 'B', 'Species', 'Drug', 'Area', 'Dataset'])
    else:
        writer.writerow(['Gene', 'baseMean', 'log2FoldChange', 'lfcSE', 'stat',
                         'pvalue', 'padj', 'Species', 'Drug', 'Area', 'Dataset'])

    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download2(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    fc = request.session.get('fc')
    print(p)
    print(fc)
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if p == '0.05' and fc == '1.2':
        sql1 = "SELECT * FROM " + str(table_name) + "_ppi"
    else:
        sql1 = "SELECT * FROM `" + str(table_name) + "_ppi_" + str(p) + "_" + str(fc) + "`"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_PPI.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['node1', 'node2', 'combined_score'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download3(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    fc = request.session.get('fc')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if p == '0.05' and fc == '1.2':
        sql1 = "SELECT * FROM " + str(table_name) + "_centrality"
    else:
        sql1 = "SELECT * FROM `" + str(table_name) + "_centrality_" + str(p) + "_" + str(fc) + "`"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_Centrality.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['node_name', 'MCC', 'MNC', 'Degree', 'EPC',
                     'BottleNeck', 'Closeness', 'Radiality', 'Betweenness',
                     'Stress', 'Katzcent', 'laplacian', 'semilocal'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download4(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    fc = request.session.get('fc')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if p == '0.05' and fc == '1.2':
        sql1 = "SELECT * FROM " + str(table_name) + "_go"
    else:
        sql1 = "SELECT * FROM `" + str(table_name) + "_go_" + str(p) + "_" + str(fc) + "`"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_Go.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['ONTOLOGY', 'ID', 'Description', 'GeneRatio', 'BgRatio',
                     'pvalue', 'p.adjust', 'qvalue', 'GeneName', 'Count'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download5(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    fc = request.session.get('fc')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if p == '0.05' and fc == '1.2':
        sql1 = "SELECT * FROM " + str(table_name) + "_kegg"
    else:
        sql1 = "SELECT * FROM `" + str(table_name) + "_kegg_" + str(p) + "_" + str(fc) + "`"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_Kegg.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['ID', 'Description', 'GeneRatio', 'BgRatio', 'pvalue',
                     'p.adjust', 'qvalue', 'GeneName', 'Count'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download6(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    deltaBeta = request.session.get('deltaBeta')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    sql1 = "SELECT * FROM methylation_" + str(table_name) + " WHERE `P.Value` < " + p + " And ABS(`deltaBeta`) >" + deltaBeta
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_DMP.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['Probe ID', 'logFC', 'AveExpr', 't', 'P.Value', 'B', 'deltaBeta',
                     'gene', 'feature', 'cgi', 'UCSC_CpG_Islands_Name', 'Drug', 'Area', 'Dataset'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download7(request):
    table_name = request.session.get('table_name')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    sql1 = "SELECT * FROM m_" + str(table_name) + "_dmr_b"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_DMR_Bumphunter.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['BumphunterDMR.seqnames', 'BumphunterDMR.start', 'BumphunterDMR.end', 'BumphunterDMR.width',
                     'BumphunterDMR.strand','BumphunterDMR.value', 'BumphunterDMR.area', 'BumphunterDMR.cluster',
                     'BumphunterDMR.indexStart', 'BumphunterDMR.indexEnd', 'BumphunterDMR.L','BumphunterDMR.clusterL',
                     'BumphunterDMR.p.value', 'BumphunterDMR.fwer', 'BumphunterDMR.p.valueArea', 'BumphunterDMR.fwerArea'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response


def download8(request):
    table_name = request.session.get('table_name')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    sql1 = "SELECT * FROM m_" + str(table_name) + "_dmr_p"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_DMR_ProbeLasso.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['ProbeLassoDMR.seqnames', 'ProbeLassoDMR.start', 'ProbeLassoDMR.end', 'ProbeLassoDMR.width',
                     'ProbeLassoDMR.strand', 'ProbeLassoDMR.dmrNo', 'ProbeLassoDMR.dmrP', 'ProbeLassoDMR.dmrpRank',
                     'ProbeLassoDMR.dmrChrom', 'ProbeLassoDMR.dmrStart', 'ProbeLassoDMR.dmrEnd', 'ProbeLassoDMR.dmrSize',
                     'ProbeLassoDMR.dmrCoreStart', 'ProbeLassoDMR.dmrCoreEnd', 'ProbeLassoDMR.dmrCoreSize',
                     'ProbeLassoDMR.ensemblID', 'ProbeLassoDMR.geneSymbol', 'ProbeLassoDMR.betaAv_Control',
                     'ProbeLassoDMR.betaAv_Nicotine'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download9(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    deltaBeta = request.session.get('deltaBeta')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if p == '0.05' and deltaBeta == '0.05':
        sql1 = "SELECT * FROM m_" + str(table_name) + "_go"
    else:
        sql1 = "SELECT * FROM `m_" + str(table_name) + "_go_" + str(p) + "_" + str(deltaBeta) + "`"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_Go.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['ID', 'Description', 'Size', 'pvalue', 'padj', 'ONTOLOGY'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def download10(request):
    table_name = request.session.get('table_name')
    p = request.session.get('p')
    deltaBeta = request.session.get('deltaBeta')
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if p == '0.05' and deltaBeta == '0.05':
        sql1 = "SELECT * FROM m_" + str(table_name) + "_kegg"
    else:
        sql1 = "SELECT * FROM `m_" + str(table_name) + "_kegg_" + str(p) + "_" + str(deltaBeta) + "`"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()

    filename = table_name + '_Kegg.csv'
    response = HttpResponse(content_type='text/csv')
    response['Content-Disposition'] = 'attachment; filename="{}"'.format(filename)
    writer = csv.writer(response)

    row = []
    writer.writerow(['ID', 'Description', 'Size', 'pvalue', 'padj'])
    for i in data1:
        for j in i:
            row.append(j)
        writer.writerow(row)
        row = []
    return response

def contact(request):

    return render(request, 'contact.html')