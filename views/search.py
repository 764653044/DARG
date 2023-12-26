from django.core import serializers
from django.http import JsonResponse, HttpResponse
from django.shortcuts import render, redirect
from django.core.paginator import Paginator
from django.views import View
from drug_addiction.models import SearchCriteria, methylationSearch
from drug_addiction import models
import pymysql
import math

config = {
    'host': '127.0.0.1',
    'port': 3306,
    'user': 'root',
    'password': '123456',
    'database': 'database',
    'charset': 'utf8',
}
def datasetSearch(request):
    """ 搜索 """
    search_data1 = request.GET.get('a')
    search_data2 = request.GET.get('b')
    search_data3 = request.GET.get('c')
    search_data4 = request.GET.get('d')
    search_data0 = request.GET.get('filter')
    p = request.GET.get('p')
    fc = request.GET.get('fc')
    if p == None:
        p = 0.05
        fc = 1.2

    if search_data1 == None:
        return render(request, 'datasetSearch.html')

    db = pymysql.connect(**config)
    cursor0 = db.cursor()
    sql0 = "SELECT type FROM dataset_type WHERE dataset = " + "'" + str(search_data4) + "'"
    cursor0.execute(sql0)
    type_dataset = cursor0.fetchall()
    type_dataset = type_dataset[0][0]
    db.close()

    if type_dataset == 'limma':
        sql11 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `P.Value` < 0.05 And ABS(`logFC`) > log2(1.2) "
        sql12 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `P.Value` < 0.05 And ABS(`logFC`) > log2(1.5) "
        sql13 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `P.Value` < 0.05 And ABS(`logFC`) > log2(2) "
        sql21 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `P.Value` < 0.01 And ABS(`logFC`) > log2(1.2) "
        sql22 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `P.Value` < 0.01 And ABS(`logFC`) > log2(1.5) "
        sql23 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `P.Value` < 0.01 And ABS(`logFC`) > log2(2) "

    else:
        sql11 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `pvalue` < 0.05 And ABS(`log2FoldChange`) > log2(1.2) "
        sql12 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `pvalue` < 0.05 And ABS(`log2FoldChange`) > log2(1.5) "
        sql13 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `pvalue` < 0.05 And ABS(`log2FoldChange`) > log2(2) "
        sql21 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `pvalue` < 0.01 And ABS(`log2FoldChange`) > log2(1.2) "
        sql22 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `pvalue` < 0.01 And ABS(`log2FoldChange`) > log2(1.5) "
        sql23 = "SELECT * FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
            search_data4) + " WHERE `pvalue` < 0.01 And ABS(`log2FoldChange`) > log2(2) "

    db = pymysql.connect(**config)
    cursor11 = db.cursor()
    cursor11.execute(sql11)
    data11 = cursor11.fetchall()
    db.close()
    n11 = len(data11)

    db = pymysql.connect(**config)
    cursor12 = db.cursor()
    cursor12.execute(sql12)
    data12 = cursor12.fetchall()
    db.close()
    n12 = len(data12)

    db = pymysql.connect(**config)
    cursor13 = db.cursor()
    cursor13.execute(sql13)
    data13 = cursor13.fetchall()
    db.close()
    n13 = len(data13)

    db = pymysql.connect(**config)
    cursor21 = db.cursor()
    cursor21.execute(sql21)
    data21 = cursor21.fetchall()
    db.close()
    n21 = len(data21)

    db = pymysql.connect(**config)
    cursor22 = db.cursor()
    cursor22.execute(sql22)
    data22 = cursor22.fetchall()
    db.close()
    n22 = len(data22)

    db = pymysql.connect(**config)
    cursor23 = db.cursor()
    cursor23.execute(sql23)
    data23 = cursor23.fetchall()
    db.close()
    n23 = len(data23)
    """ 获取火山图上调基因 """
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if type_dataset == 'limma':
        if p == '0.01' and fc == '1.2':
            sql1 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > log2(1.2) AND `P.Value` <0.01 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.01' and fc == '1.5':
            sql1 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > log2(1.5) AND `P.Value` <0.01 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.05' and fc == '1.5':
            sql1 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > log2(1.5) AND `P.Value` <0.05 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.01' and fc == '2':
            sql1 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > log2(2) AND `P.Value` <0.01 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.05' and fc == '2':
            sql1 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > log2(2) AND `P.Value` <0.05 AND `P.Value` != 'NA' AND `P.Value` != 0"
        else:
            sql1 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > log2(1.2) AND `P.Value` <0.05 AND `P.Value` != 'NA' AND `P.Value` != 0"
    else:
        if p == '0.01' and fc == '1.2':
            sql1 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > log2(1.2) AND `pvalue` <0.01 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.01' and fc == '1.5':
            sql1 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > log2(1.5) AND `pvalue` <0.01 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.05' and fc == '1.5':
            sql1 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > log2(1.5) AND `pvalue` <0.05 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.01' and fc == '2':
            sql1 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > log2(2) AND `pvalue` <0.01 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.05' and fc == '2':
            sql1 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > log2(2) AND `pvalue` <0.05 AND `pvalue` != 'NA' AND `pvalue` != 0"
        else:
            sql1 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > log2(1.2) AND `pvalue` <0.05 AND `pvalue` != 'NA' AND `pvalue` != 0"

    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()
    table_list1 = []
    z_list1 = []
    for z1 in data1:
        for z2 in z1:
            z_list1.append(z2)
            float_list1 = [float(x) for x in z_list1]
        table_list1.append(float_list1)
        z_list1 = []
    for i in table_list1:
        i[1] = -math.log10(i[1])
    """ 获取火山图下调基因 """
    db = pymysql.connect(**config)
    cursor2 = db.cursor()
    if type_dataset == 'limma':
        if p == '0.01' and fc == '1.2':
            sql2 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` < -log2(1.2) AND `P.Value` <0.01 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.01' and fc == '1.5':
            sql2 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` < -log2(1.5) AND `P.Value` <0.01 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.05' and fc == '1.5':
            sql2 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` < -log2(1.5) AND `P.Value` <0.05 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.01' and fc == '2':
            sql2 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` < -log2(2) AND `P.Value` <0.01 AND `P.Value` != 'NA' AND `P.Value` != 0"
        elif p == '0.05' and fc == '2':
            sql2 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` < -log2(2) AND `P.Value` <0.05 AND `P.Value` != 'NA' AND `P.Value` != 0"
        else:
            sql2 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` < -log2(1.2) AND `P.Value` <0.05 AND `P.Value` != 'NA' AND `P.Value` != 0"
    else:
        if p == '0.01' and fc == '1.2':
            sql2 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` < -log2(1.2) AND `pvalue` <0.01 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.01' and fc == '1.5':
            sql2 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` < -log2(1.5) AND `pvalue` <0.01 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.05' and fc == '1.5':
            sql2 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` < -log2(1.5) AND `pvalue` <0.05 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.01' and fc == '2':
            sql2 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` < -log2(2) AND `pvalue` <0.01 AND `pvalue` != 'NA' AND `pvalue` != 0"
        elif p == '0.05' and fc == '2':
            sql2 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` < -log2(2) AND `pvalue` <0.05 AND `pvalue` != 'NA' AND `pvalue` != 0"
        else:
            sql2 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` < -log2(1.2) AND `pvalue` <0.05 AND `pvalue` != 'NA' AND `pvalue` != 0"
    cursor2.execute(sql2)
    data2 = cursor2.fetchall()
    db.close()
    table_list2 = []
    z_list2 = []
    for z1 in data2:
        for z2 in z1:
            z_list2.append(z2)
            float_list2 = [float(x) for x in z_list2]
        table_list2.append(float_list2)
        z_list2 = []
    for i in table_list2:
        i[1] = -math.log10(i[1])

    """ 获取火山图不差异基因 """
    db = pymysql.connect(**config)
    cursor3 = db.cursor()
    if type_dataset == 'limma':
        if p == '0.01' and fc == '1.2':
            sql3 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > -log2(1.2) AND `logFC` < log2(1.2) AND `P.Value` != 'NA' AND `P.Value` != 0 OR `P.Value` >0.01"
        elif p == '0.01' and fc == '1.5':
            sql3 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > -log2(1.5) AND `logFC` < log2(1.5) AND `P.Value` != 'NA' AND `P.Value` != 0 OR `P.Value` >0.01"
        elif p == '0.05' and fc == '1.5':
            sql3 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > -log2(1.5) AND `logFC` < log2(1.5) AND `P.Value` != 'NA' AND `P.Value` != 0 OR `P.Value` >0.05"
        elif p == '0.01' and fc == '2':
            sql3 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > -log2(2) AND `logFC` < log2(2) AND `P.Value` != 'NA' AND `P.Value` != 0 OR `P.Value` >0.01"
        elif p == '0.05' and fc == '2':
            sql3 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > -log2(2) AND `logFC` < log2(2) AND `P.Value` != 'NA' AND `P.Value` != 0 OR `P.Value` >0.05"
        else:
            sql3 = "SELECT `logFC`,`P.Value` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `logFC` > -log2(1.2) AND `logFC` < log2(1.2) AND `P.Value` != 'NA' AND `P.Value` != 0 OR `P.Value` >0.05"
    else:
        if p == '0.01' and fc == '1.2':
            sql3 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > -log2(1.2) AND `log2FoldChange` < log2(1.2) AND `pvalue` != 'NA' AND `pvalue` != 0 OR `pvalue` >0.01"
        elif p == '0.01' and fc == '1.5':
            sql3 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > -log2(1.5) AND `log2FoldChange` < log2(1.5) AND `pvalue` != 'NA' AND `pvalue` != 0 OR `pvalue` >0.01"
        elif p == '0.05' and fc == '1.5':
            sql3 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > -log2(1.5) AND `log2FoldChange` < log2(1.5) AND `pvalue` != 'NA' AND `pvalue` != 0 OR `pvalue` >0.05"
        elif p == '0.01' and fc == '2':
            sql3 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > -log2(2) AND `log2FoldChange` < log2(2) AND `pvalue` != 'NA' AND `pvalue` != 0 OR `pvalue` >0.01"
        elif p == '0.05' and fc == '2':
            sql3 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > -log2(2) AND `log2FoldChange` < log2(2) AND `pvalue` != 'NA' AND `pvalue` != 0 OR `pvalue` >0.05"
        else:
            sql3 = "SELECT `log2FoldChange`,`pvalue` FROM " + str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `log2FoldChange` > -log2(1.2) AND `log2FoldChange` < log2(1.2) AND `pvalue` != 'NA' AND `pvalue` != 0 OR `pvalue` >0.05"
    cursor3.execute(sql3)
    data3 = cursor3.fetchall()
    db.close()
    table_list3 = []
    z_list3 = []
    for z1 in data3:
        for z2 in z1:
            z_list3.append(z2)
            float_list3 = [float(x) for x in z_list3]
        table_list3.append(float_list3)
        z_list3 = []
    for i in table_list3:
        i[1] = -math.log10(i[1])

    """ 获取PPI数据 """
    db = pymysql.connect(**config)
    cursor4 = db.cursor()
    if p == '0.01' and fc == '1.2':
        sql4 = "SELECT `protein1`,`protein2` FROM `" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_ppi_0.01_1.2`" + " WHERE `combined_score` > 700 OR `combined_score` = 700"
    elif p == '0.01' and fc == '1.5':
        sql4 = "SELECT `protein1`,`protein2` FROM `" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_ppi_0.01_1.5`" + " WHERE `combined_score` > 700 OR `combined_score` = 700"
    elif p == '0.05' and fc == '1.5':
        sql4 = "SELECT `protein1`,`protein2` FROM `" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_ppi_0.05_1.5`" + " WHERE `combined_score` > 700 OR `combined_score` = 700"
    elif p == '0.01' and fc == '2':
        sql4 = "SELECT `protein1`,`protein2` FROM `" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_ppi_0.01_2`" + " WHERE `combined_score` > 700 OR `combined_score` = 700"
    elif p == '0.05' and fc == '2':
        sql4 = "SELECT `protein1`,`protein2` FROM `" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_ppi_0.05_2`" + " WHERE `combined_score` > 700 OR `combined_score` = 700"
    else:
        sql4 = "SELECT `protein1`,`protein2` FROM `" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_ppi`" + " WHERE `combined_score` > 700 OR `combined_score` = 700"
    cursor4.execute(sql4)
    ppi = cursor4.fetchall()
    db.close()
    ppi_data1 = []
    ppi_data2 = []
    for j in ppi:
        ppi_data1.append(j[0])
        ppi_data2.append(j[1])
    zipped_ppi = zip(ppi_data1, ppi_data2)

    """ 获取中心性算法得分数据 """
    if p == '0.01' and fc == '1.2':
        centrality_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_centrality_0.01_1.2"
        go_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_go_0.01_1.2"
        kegg_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_kegg_0.01_1.2"
    elif p == '0.01' and fc == '1.5':
        centrality_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_centrality_0.01_1.5"
        go_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_go_0.01_1.5"
        kegg_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_kegg_0.01_1.5"
    elif p == '0.05' and fc == '1.5':
        centrality_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_centrality_0.05_1.5"
        go_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_go_0.05_1.5"
        kegg_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_kegg_0.05_1.5"
    elif p == '0.01' and fc == '2':
        centrality_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_centrality_0.01_2"
        go_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_go_0.01_2"
        kegg_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_kegg_0.01_2"
    elif p == '0.05' and fc == '2':
        centrality_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_centrality_0.05_2"
        go_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_go_0.05_2"
        kegg_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_kegg_0.05_2"
    else:
        centrality_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_centrality"
        go_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_go"
        kegg_table_name = str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_" + str(search_data4) + "_kegg"

    db = pymysql.connect(**config)
    cursor24 = db.cursor()
    sql24 = "SELECT `node_name`,`Degree` FROM `" + str(centrality_table_name) + "`"
    cursor24.execute(sql24)
    degree_all = cursor24.fetchall()
    db.close()
    centrality_length = len(degree_all)
    centrality_number = centrality_length//10
    # if centrality_length > 200:
    #     centrality_number = 30
    # elif centrality_length < 200 and centrality_length > 100:
    #     centrality_number = 20
    # elif centrality_length < 100 and centrality_length > 50:
    #     centrality_number = 10
    # else:
    #     centrality_number = 5
    """ 获取中心性算法得分数据 """
    db = pymysql.connect(**config)
    cursor5 = db.cursor()
    sql5 = "SELECT `node_name`,`Degree` FROM `" + str(centrality_table_name) + "` ORDER BY `Degree` DESC LIMIT " + str(centrality_number)
    cursor5.execute(sql5)
    degree = cursor5.fetchall()
    db.close()
    degree_list = []
    degree_list1 = []
    degree_name = []
    for i in degree:
        degree_list1.append(i[0])
        degree_list1.append(i[1])
        degree_list.append(degree_list1)
        degree_list1 = []
        degree_name.append(i[0])
    set_degree = set(degree_name)

    db = pymysql.connect(**config)
    cursor6 = db.cursor()
    sql6 = "SELECT `node_name`,`Radiality` FROM `" + str(centrality_table_name) + "` ORDER BY `Radiality` DESC LIMIT " + str(centrality_number)
    cursor6.execute(sql6)
    Radiality = cursor6.fetchall()
    db.close()
    Radiality_list = []
    Radiality_list1 = []
    Radiality_name = []
    for i in Radiality:
        Radiality_list1.append(i[0])
        Radiality_list1.append(i[1])
        Radiality_list.append(Radiality_list1)
        Radiality_list1 = []
        Radiality_name.append(i[0])
    set_Radiality = set(Radiality_name)

    # db = pymysql.connect(**config)
    # cursor7 = db.cursor()
    # sql7 = "SELECT `node_name`,`Closeness` FROM `" + str(centrality_table_name) + "` ORDER BY `Closeness` DESC LIMIT " + str(centrality_number)
    # cursor7.execute(sql7)
    # Closeness = cursor7.fetchall()
    # db.close()
    # Closeness_list = []
    # Closeness_list1 = []
    # Closeness_name = []
    # for i in Closeness:
    #     Closeness_list1.append(i[0])
    #     Closeness_list1.append(i[1])
    #     Closeness_list.append(Closeness_list1)
    #     Closeness_list1 = []
    #     Closeness_name.append(i[0])
    # set_Closeness = set(Closeness_name)
    #
    # db = pymysql.connect(**config)
    # cursor8 = db.cursor()
    # sql8 = "SELECT `node_name`,`MCC` FROM `" + str(centrality_table_name) + "` ORDER BY `MCC` DESC LIMIT " + str(centrality_number)
    # cursor8.execute(sql8)
    # MCC = cursor8.fetchall()
    # db.close()
    # MCC_list = []
    # MCC_list1 = []
    # MCC_name = []
    # for i in Closeness:
    #     MCC_list1.append(i[0])
    #     MCC_list1.append(i[1])
    #     MCC_list.append(MCC_list1)
    #     MCC_list1 = []
    #     MCC_name.append(i[0])
    # set_MCC = set(MCC_name)

    # db = pymysql.connect(**config)
    # cursor9 = db.cursor()
    # sql9 = "SELECT `node_name`,`Betweenness` FROM `" + str(centrality_table_name) + "` ORDER BY `Betweenness` DESC LIMIT " + str(centrality_number)
    # cursor9.execute(sql9)
    # Betweenness = cursor9.fetchall()
    # db.close()
    # Betweenness_list = []
    # Betweenness_list1 = []
    # Betweenness_name = []
    # for i in Betweenness:
    #     Betweenness_list1.append(i[0])
    #     Betweenness_list1.append(i[1])
    #     Betweenness_list.append(Betweenness_list1)
    #     Betweenness_list1 = []
    #     Betweenness_name.append(i[0])
    # set_Betweenness = set(Betweenness_name)

    db = pymysql.connect(**config)
    cursor10 = db.cursor()
    sql10 = "SELECT `node_name`,`EPC` FROM `" + str(centrality_table_name) + "` ORDER BY `EPC` DESC LIMIT " + str(centrality_number)
    cursor10.execute(sql10)
    EPC = cursor10.fetchall()
    db.close()
    EPC_list = []
    EPC_list1 = []
    EPC_name = []
    for i in EPC:
        EPC_list1.append(i[0])
        EPC_list1.append(i[1])
        EPC_list.append(EPC_list1)
        EPC_list1 = []
        EPC_name.append(i[0])
    set_EPC = set(EPC_name)

    db = pymysql.connect(**config)
    cursor14 = db.cursor()
    sql14 = "SELECT `node_name`,`MNC` FROM `" + str(centrality_table_name) + "` ORDER BY `MNC` DESC LIMIT " + str(centrality_number)
    cursor14.execute(sql14)
    MNC = cursor14.fetchall()
    db.close()
    MNC_list = []
    MNC_list1 = []
    MNC_name = []
    for i in MNC:
        MNC_list1.append(i[0])
        MNC_list1.append(i[1])
        MNC_list.append(MNC_list1)
        MNC_list1 = []
        MNC_name.append(i[0])
    set_MNC = set(MNC_name)

    db = pymysql.connect(**config)
    cursor16 = db.cursor()
    sql16 = "SELECT `node_name`,`Katzcent` FROM `" + str(centrality_table_name) + "` ORDER BY `Katzcent` DESC LIMIT " + str(centrality_number)
    cursor16.execute(sql16)
    Katz = cursor16.fetchall()
    db.close()
    Katz_list = []
    Katz_list1 = []
    Katz_name = []
    for i in Katz:
        Katz_list1.append(i[0])
        Katz_list1.append(i[1])
        Katz_list.append(Katz_list1)
        Katz_list1 = []
        Katz_name.append(i[0])
    set_Katz = set(Katz_name)

    # db = pymysql.connect(**config)
    # cursor15 = db.cursor()
    # sql15 = "SELECT `node_name`,`Stress` FROM `" + str(centrality_table_name) + "` ORDER BY `Stress` DESC LIMIT " + str(centrality_number)
    # cursor15.execute(sql15)
    # Stress = cursor15.fetchall()
    # db.close()
    # Stress_list = []
    # Stress_list1 = []
    # Stress_name = []
    # for i in Stress:
    #     Stress_list1.append(i[0])
    #     Stress_list1.append(i[1])
    #     Stress_list.append(Stress_list1)
    #     Stress_list1 = []
    #     Stress_name.append(i[0])
    # set_Stress = set(Stress_name)

    db = pymysql.connect(**config)
    cursor17 = db.cursor()
    sql17 = "SELECT `node_name`,`laplacian` FROM `" + str(centrality_table_name) + "` ORDER BY `laplacian` DESC LIMIT " + str(centrality_number)
    cursor17.execute(sql17)
    Laplacian = cursor17.fetchall()
    db.close()
    Laplacian_list = []
    Laplacian_list1 = []
    Laplacian_name = []
    for i in Laplacian:
        Laplacian_list1.append(i[0])
        Laplacian_list1.append(i[1])
        Laplacian_list.append(Laplacian_list1)
        Laplacian_list1 = []
        Laplacian_name.append(i[0])
    set_Laplacian = set(Laplacian_name)

    db = pymysql.connect(**config)
    cursor18 = db.cursor()
    sql18 = "SELECT `node_name`,`semilocal` FROM `" + str(centrality_table_name) + "` ORDER BY `semilocal` DESC LIMIT " + str(centrality_number)
    cursor18.execute(sql18)
    SLC = cursor18.fetchall()
    db.close()
    SLC_list = []
    SLC_list1 = []
    SLC_name = []
    for i in SLC:
        SLC_list1.append(i[0])
        SLC_list1.append(i[1])
        SLC_list.append(SLC_list1)
        SLC_list1 = []
        SLC_name.append(i[0])
    set_SLC = set(SLC_name)

    # db = pymysql.connect(**config)
    # cursor19 = db.cursor()
    # sql19 = "SELECT `node_name`,`BottleNeck` FROM `" + str(centrality_table_name) + "` ORDER BY `BottleNeck` DESC LIMIT " + str(centrality_number)
    # cursor19.execute(sql19)
    # BottleNeck  = cursor19.fetchall()
    # db.close()
    # BottleNeck_list = []
    # BottleNeck_list1 = []
    # BottleNeck_name = []
    # for i in BottleNeck:
    #     BottleNeck_list1.append(i[0])
    #     BottleNeck_list1.append(i[1])
    #     BottleNeck_list.append(BottleNeck_list1)
    #     BottleNeck_list1 = []
    #     BottleNeck_name.append(i[0])
    # set_BottleNeck = set(BottleNeck_name)

    # intersection = set_degree & set_Radiality & set_Closeness & set_MCC & set_Betweenness & set_EPC & set_MNC & set_Katz & set_Stress & set_Laplacian & set_SLC & set_BottleNeck
    intersection = set_degree & set_Radiality & set_EPC & set_MNC & set_Katz & set_Laplacian & set_SLC
    intersection_list = []
    for i in intersection:
        intersection_list.append(i)
    if centrality_length < 30:
        intersection_list = []
    request.session['intersection_list'] = intersection_list

    # go-BP富集分析
    db = pymysql.connect(**config)
    cursor21 = db.cursor()
    sql21 = "SELECT `Count`,`Description`,`pvalue`,`p.adjust`,`qvalue`,`geneName`,`ID` FROM `" + go_table_name + "` WHERE ONTOLOGY = 'BP' LIMIT 10"
    cursor21.execute(sql21)
    go_bp = cursor21.fetchall()
    db.close()
    BP_list = []
    BP_list1 = []
    GO_name = []
    GO_count = []
    GO_padjust = []
    GO_list = []
    for z in go_bp:
        GO_count.append(z[0])
        GO_name.append(z[1])
        GO_padjust.append(z[3])
        for y in z:
            BP_list1.append(y)
        BP_list.append(BP_list1)
        GO_list.append(BP_list1)
        BP_list1 = []

    # go-CC富集分析
    db = pymysql.connect(**config)
    cursor22 = db.cursor()
    sql22 = "SELECT `Count`,`Description`,`pvalue`,`p.adjust`,`qvalue`,`geneName`,`ID` FROM `" + go_table_name + "` WHERE ONTOLOGY = 'CC' LIMIT 10"
    cursor22.execute(sql22)
    go_cc = cursor22.fetchall()
    db.close()
    CC_list = []
    CC_list1 = []
    for z in go_cc:
        GO_count.append(z[0])
        GO_name.append(z[1])
        GO_padjust.append(z[3])
        for y in z:
            CC_list1.append(y)
        CC_list.append(CC_list1)
        GO_list.append(CC_list1)
        CC_list1 = []

    # go-MF富集分析
    db = pymysql.connect(**config)
    cursor23 = db.cursor()
    sql23 = "SELECT `Count`,`Description`,`pvalue`,`p.adjust`,`qvalue`,`geneName`,`ID` FROM `" + go_table_name + "` WHERE ONTOLOGY = 'MF' LIMIT 10"
    cursor23.execute(sql23)
    go_mf = cursor23.fetchall()
    db.close()
    MF_list = []
    MF_list1 = []
    for z in go_mf:
        GO_count.append(z[0])
        GO_name.append(z[1])
        GO_padjust.append(z[3])
        for y in z:
            MF_list1.append(y)
        MF_list.append(MF_list1)
        GO_list.append(MF_list1)
        MF_list1 = []

    if GO_count != []:
        GO_count_max = max(GO_count)
        GO_count_min = min(GO_count)
        GO_padjust_max = max(GO_padjust)
        GO_padjust_min = min(GO_padjust)
        GO_name.reverse()

    # kegg富集分析
    db = pymysql.connect(**config)
    cursor20 = db.cursor()
    sql20 = "SELECT `Count`,`Description`,`pvalue`,`p.adjust`,`qvalue`,`geneName`,`ID` FROM `" + kegg_table_name + "` LIMIT 30"
    cursor20.execute(sql20)
    kegg = cursor20.fetchall()
    db.close()
    kegg_list = []
    kegg_list1 = []
    kegg_name = []
    kegg_count = []
    kegg_padjust = []
    for z in kegg:
        kegg_count.append(z[0])
        kegg_name.append(z[1])
        kegg_padjust.append(z[3])
        for y in z:
            kegg_list1.append(y)
        kegg_list.append(kegg_list1)
        kegg_list1 = []
    if kegg != ():
        kegg_count_max = max(kegg_count)
        kegg_count_min = min(kegg_count)
        kegg_padjust_max = max(kegg_padjust)
        kegg_padjust_min = min(kegg_padjust)
        kegg_name.reverse()

    return render(request, 'datasetSearch.html', locals())

def degs(request):
    """ 搜索 """
    search_data1 = request.GET.get('a')
    search_data2 = request.GET.get('b')
    search_data3 = request.GET.get('c')
    search_data4 = request.GET.get('d')
    p = request.GET.get('p')
    fc = request.GET.get('fc')
    search_data0 = request.GET.get('filter')
    key = request.GET.get('key')
    key_genes = request.session.get('intersection_list')
    sort = request.GET.get('sort')
    dir = request.GET.get('dir')

    db = pymysql.connect(**config)
    cursor2 = db.cursor()
    sql2 = "SELECT type FROM dataset_type WHERE dataset = " + "'" + str(search_data4) + "'"
    cursor2.execute(sql2)
    type_dataset = cursor2.fetchall()
    type_dataset = type_dataset[0][0]
    db.close()

    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if search_data0 != '' and search_data0 != None:
        if type_dataset == 'limma':
            if sort == 'Pvalue':
                sotr = 'P.Value'
            elif sort == 'adjPvalue':
                sotr = 'adj.P.Val'
            else:
                sotr = 'logFC'

            if dir == 'desc':
                sql = "SELECT `Gene`,ROUND(`logFC`, 8) AS `logFC`,ROUND(`AveExpr`, 8) AS `AveExpr`,ROUND(`t`, 8) AS `t`," \
                      "ROUND(`P.Value`, 8) AS `P.Value`,ROUND(`adj.P.Val`, 8) AS `adj.P.Val`,ROUND(`B`, 8) AS `B` FROM " + \
                      str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ") And Gene LIKE " + "'%" + str(
                    search_data0) + "%' ORDER BY CONVERT(`" + str(sotr) +"`, DOUBLE) DESC;"
            elif dir == 'asc':
                sql = "SELECT `Gene`,ROUND(`logFC`, 8) AS `logFC`,ROUND(`AveExpr`, 8) AS `AveExpr`,ROUND(`t`, 8) AS `t`," \
                      "ROUND(`P.Value`, 8) AS `P.Value`,ROUND(`adj.P.Val`, 8) AS `adj.P.Val`,ROUND(`B`, 8) AS `B` FROM " + \
                      str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ") And Gene LIKE " + "'%" + str(
                    search_data0) + "%' ORDER BY CONVERT(`" + str(sotr) +"`, DOUBLE) ASC;"
            else:
                sql = "SELECT `Gene`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                      "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                      "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                      "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                      "CASE WHEN ABS(`adj.P.Val`) < 1E-8 THEN CAST(`adj.P.Val` AS CHAR) ELSE CAST(`adj.P.Val` AS DECIMAL(18, 8)) END AS `adj.P.Val`," \
                      "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B` FROM " + \
                      str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ") And Gene LIKE " + "'%" + str(
                    search_data0) + "%'"
        else:
            if sort == 'Pvalue':
                sotr = 'pvalue'
            elif sort == 'adjPvalue':
                sotr = 'padj'
            else:
                sotr = 'log2FoldChange'

            if dir == 'desc':
                sql = "SELECT `Gene`,ROUND(`baseMean`, 8) AS `baseMean`,ROUND(`log2FoldChange`, 8) AS `log2FoldChange`,ROUND(`lfcSE`, 8) " \
                      "AS `lfcSE`,ROUND(`stat`, 8) AS `stat`,ROUND(`pvalue`, 8) AS `pvalue`,ROUND(`padj`, 8) AS `padj` FROM " + str(
                    search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ") And Gene LIKE " + "'%" + str(
                    search_data0) + "%' ORDER BY CONVERT(`" + str(sotr) +"`, DOUBLE) DESC;"
            elif dir == 'asc':
                sql = "SELECT `Gene`,ROUND(`baseMean`, 8) AS `baseMean`,ROUND(`log2FoldChange`, 8) AS `log2FoldChange`,ROUND(`lfcSE`, 8) " \
                      "AS `lfcSE`,ROUND(`stat`, 8) AS `stat`,ROUND(`pvalue`, 8) AS `pvalue`,ROUND(`padj`, 8) AS `padj` FROM " + str(
                    search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ") And Gene LIKE " + "'%" + str(
                    search_data0) + "%' ORDER BY CONVERT(`" + str(sotr) +"`, DOUBLE) ASC;"
            else:
                sql = "SELECT `Gene`,CASE WHEN ABS(`baseMean`) < 1E-8 THEN CAST(`baseMean` AS CHAR) ELSE CAST(`baseMean` AS DECIMAL(18, 8)) END AS `baseMean`," \
                      "CASE WHEN ABS(`log2FoldChange`) < 1E-8 THEN CAST(`log2FoldChange` AS CHAR) ELSE CAST(`log2FoldChange` AS DECIMAL(18, 8)) END AS `log2FoldChange`," \
                      "CASE WHEN ABS(`lfcSE`) < 1E-8 THEN CAST(`lfcSE` AS CHAR) ELSE CAST(`lfcSE` AS DECIMAL(18, 8)) END AS `lfcSE`," \
                      "CASE WHEN ABS(`stat`) < 1E-8 THEN CAST(`stat` AS CHAR) ELSE CAST(`stat` AS DECIMAL(18, 8)) END AS `stat`," \
                      "CASE WHEN ABS(`pvalue`) < 1E-8 THEN CAST(`pvalue` AS CHAR) ELSE CAST(`pvalue` AS DECIMAL(18, 8)) END AS `pvalue`," \
                      "CASE WHEN ABS(`padj`) < 1E-8 THEN CAST(`padj` AS CHAR) ELSE CAST(`padj` AS DECIMAL(18, 8)) END AS `padj` FROM " + str(
                    search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ") And Gene LIKE " + "'%" + str(
                    search_data0) + "%'"
    else:
        if type_dataset == 'limma':
            if sort == 'Pvalue':
                sotr = 'P.Value'
            elif sort == 'adjPvalue':
                sotr = 'adj.P.Val'
            else:
                sotr = 'logFC'

            if dir == 'desc':
                sql = "SELECT `Gene`,ROUND(`logFC`, 8) AS `logFC`,ROUND(`AveExpr`, 8) AS `AveExpr`,ROUND(`t`, 8) AS `t`," \
                      "ROUND(`P.Value`, 8) AS `P.Value`,ROUND(`adj.P.Val`, 8) AS `adj.P.Val`,ROUND(`B`, 8) AS `B` FROM " + \
                      str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ") ORDER BY CONVERT(`"+ str(sotr) +"`, DOUBLE) DESC;"
            elif dir == 'asc':
                sql = "SELECT `Gene`,ROUND(`logFC`, 8) AS `logFC`,ROUND(`AveExpr`, 8) AS `AveExpr`,ROUND(`t`, 8) AS `t`," \
                      "ROUND(`P.Value`, 8) AS `P.Value`,ROUND(`adj.P.Val`, 8) AS `adj.P.Val`,ROUND(`B`, 8) AS `B` FROM " + \
                      str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ") ORDER BY CONVERT(`"+ str(sotr) +"`, DOUBLE) ASC;"
            else:
                sql = "SELECT `Gene`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                      "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                      "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                      "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                      "CASE WHEN ABS(`adj.P.Val`) < 1E-8 THEN CAST(`adj.P.Val` AS CHAR) ELSE CAST(`adj.P.Val` AS DECIMAL(18, 8)) END AS `adj.P.Val`," \
                      "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B` FROM " + \
                      str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ") "
        else:
            if sort == 'Pvalue':
                sotr = 'pvalue'
            elif sort == 'adjPvalue':
                sotr = 'padj'
            else:
                sotr = 'log2FoldChange'

            if dir == 'desc':
                sql = "SELECT `Gene`,ROUND(`baseMean`, 8) AS `baseMean`,ROUND(`log2FoldChange`, 8) AS `log2FoldChange`,ROUND(`lfcSE`, 8) " \
                      "AS `lfcSE`,ROUND(`stat`, 8) AS `stat`,ROUND(`pvalue`, 8) AS `pvalue`,ROUND(`padj`, 8) AS `padj` FROM " + str(
                    search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ") ORDER BY CONVERT(`"+ str(sotr) +"`, DOUBLE) DESC;"
            elif dir == 'asc':
                sql = "SELECT `Gene`,ROUND(`baseMean`, 8) AS `baseMean`,ROUND(`log2FoldChange`, 8) AS `log2FoldChange`,ROUND(`lfcSE`, 8) " \
                      "AS `lfcSE`,ROUND(`stat`, 8) AS `stat`,ROUND(`pvalue`, 8) AS `pvalue`,ROUND(`padj`, 8) AS `padj` FROM " + str(
                    search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ") ORDER BY CONVERT(`"+ str(sotr) +"`, DOUBLE) ASC;"
            else:
                sql = "SELECT `Gene`,CASE WHEN ABS(`baseMean`) < 1E-8 THEN CAST(`baseMean` AS CHAR) ELSE CAST(`baseMean` AS DECIMAL(18, 8)) END AS `baseMean`," \
                      "CASE WHEN ABS(`log2FoldChange`) < 1E-8 THEN CAST(`log2FoldChange` AS CHAR) ELSE CAST(`log2FoldChange` AS DECIMAL(18, 8)) END AS `log2FoldChange`," \
                      "CASE WHEN ABS(`lfcSE`) < 1E-8 THEN CAST(`lfcSE` AS CHAR) ELSE CAST(`lfcSE` AS DECIMAL(18, 8)) END AS `lfcSE`," \
                      "CASE WHEN ABS(`stat`) < 1E-8 THEN CAST(`stat` AS CHAR) ELSE CAST(`stat` AS DECIMAL(18, 8)) END AS `stat`," \
                      "CASE WHEN ABS(`pvalue`) < 1E-8 THEN CAST(`pvalue` AS CHAR) ELSE CAST(`pvalue` AS DECIMAL(18, 8)) END AS `pvalue`," \
                      "CASE WHEN ABS(`padj`) < 1E-8 THEN CAST(`padj` AS CHAR) ELSE CAST(`padj` AS DECIMAL(18, 8)) END AS `padj` FROM " + str(
                    search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ") "
    cursor1.execute(sql)
    data = cursor1.fetchall()
    db.close()

    key_list = []
    for i in key_genes:
        db = pymysql.connect(**config)
        cursor4 = db.cursor()
        if type_dataset == 'limma':
            sql4 = "SELECT `Gene`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                      "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                      "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                      "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                      "CASE WHEN ABS(`adj.P.Val`) < 1E-8 THEN CAST(`adj.P.Val` AS CHAR) ELSE CAST(`adj.P.Val` AS DECIMAL(18, 8)) END AS `adj.P.Val`," \
                      "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B` FROM " + \
                      str(search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                    search_data4) + " WHERE `P.Value` < " + p + " And ABS(`logFC`) > log2(" + fc + ") And Gene = " + "'" + str(
                i) + "'"
        else:
            sql4 = "SELECT `Gene`,CASE WHEN ABS(`baseMean`) < 1E-8 THEN CAST(`baseMean` AS CHAR) ELSE CAST(`baseMean` AS DECIMAL(18, 8)) END AS `baseMean`," \
                      "CASE WHEN ABS(`log2FoldChange`) < 1E-8 THEN CAST(`log2FoldChange` AS CHAR) ELSE CAST(`log2FoldChange` AS DECIMAL(18, 8)) END AS `log2FoldChange`," \
                      "CASE WHEN ABS(`lfcSE`) < 1E-8 THEN CAST(`lfcSE` AS CHAR) ELSE CAST(`lfcSE` AS DECIMAL(18, 8)) END AS `lfcSE`," \
                      "CASE WHEN ABS(`stat`) < 1E-8 THEN CAST(`stat` AS CHAR) ELSE CAST(`stat` AS DECIMAL(18, 8)) END AS `stat`," \
                      "CASE WHEN ABS(`pvalue`) < 1E-8 THEN CAST(`pvalue` AS CHAR) ELSE CAST(`pvalue` AS DECIMAL(18, 8)) END AS `pvalue`," \
                      "CASE WHEN ABS(`padj`) < 1E-8 THEN CAST(`padj` AS CHAR) ELSE CAST(`padj` AS DECIMAL(18, 8)) END AS `padj` FROM " + str(
                    search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(
                search_data4) + " WHERE `pvalue` < " + p + " And ABS(`log2FoldChange`) > log2(" + fc + ") And Gene = " + "'" + str(
                i) + "'"
        cursor4.execute(sql4)
        key_content = cursor4.fetchall()
        db.close()
        if key_content:
            key_list.append(key_content[0])

    if type_dataset == 'limma':
        if sort == 'Pvalue':
            n = 4
        elif sort == 'adjPvalue':
            n = 5
        else:
            n = 1
    else:
        if sort == 'Pvalue':
            n = 5
        elif sort == 'adjPvalue':
            n = 6
        else:
            n = 2

    if dir == 'desc':
        sorted_data = sorted(key_list, key=lambda x: x[n], reverse=True)
    elif dir == 'asc':
        sorted_data = sorted(key_list, key=lambda x: x[n], reverse=False)
    else:
        sorted_data = key_list

    db = pymysql.connect(**config)
    cursor3 = db.cursor()
    sql = "SELECT strain FROM strain_mouse WHERE dataset = " + "'" + str(search_data4) + "'"
    cursor3.execute(sql)
    strain = cursor3.fetchall()
    db.close()
    if strain:
        strain = strain[0][0]

    if key == '1':
        data = sorted_data

    # 生成分页器对象
    paginator = Paginator(data, 10)
    # 获取浏览器端请求的页码,需要设置默认值
    current_page_num = int(request.GET.get('page', 1))
    # 当前页的所有书对象
    current_page = paginator.page(current_page_num)
    # 页码列表，可迭代
    # 分页过多，需要用条件判断显示的页码
    if paginator.num_pages > 11:  # 11就是显示的固定个数
        if current_page_num - 5 < 1:  # 接近最小页码时，固定页码，否则会出现负数
            page_range = range(1, 12)
        elif current_page_num + 5 > paginator.num_pages:  # 接近最大页码时，根据最大页码限制页码数，否则会出现不存在的页码
            page_range = range(paginator.num_pages - 10, paginator.num_pages + 1)
        else:
            page_range = range(current_page_num - 5, current_page_num + 6)  # 其他情况，显示当前的挨着的几个
    else:
        page_range = paginator.page_range  # 页码总数不足时，显示全部，即不会超宽
    return render(request, 'datasetSearch_DEGs.html', locals())


class AreaView(View):
    def get(self, request):
        return render(request, 'datasetSearch.html')

def getInfo(request):
    pid = request.GET.get('pid', -1)
    pid = int(pid)
    areaList = SearchCriteria.objects.filter(parentid=pid)
    jareaList = serializers.serialize('json', areaList)
    return JsonResponse({'jareaList': jareaList})

def getMethylationInfo(request):
    pid = request.GET.get('pid', -1)
    pid = int(pid)
    areaList = methylationSearch.objects.filter(parentid=pid)
    jareaList = serializers.serialize('json', areaList)
    return JsonResponse({'jareaList': jareaList})

def gene(request):
    gene = request.GET.get('gene')
    a = request.GET.get('a')
    b = request.GET.get('b')
    c = request.GET.get('c')
    d = request.GET.get('d')
    # 将小写变成大写
    d_upper = d.upper()
    # 获取Control值
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    sql1 = "SELECT Control FROM " + str(a) + "_" + str(b) + "_" + str(c) + "_" + str(d) + "_gene" + " WHERE Gene = " + "'" + str(gene) + "'"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    data1 = data1[0][0]
    db.close()
    # 获取Drug值
    db = pymysql.connect(**config)
    cursor2 = db.cursor()
    sql2 = "SELECT Drug FROM " + str(a) + "_" + str(b) + "_" + str(c) + "_" + str(d) + "_gene" + " WHERE Gene = " + "'" + str(gene) + "'"
    cursor2.execute(sql2)
    data2 = cursor2.fetchall()
    data2 = data2[0][0]
    db.close()
    # 获取横坐标值
    db = pymysql.connect(**config)
    cursor3 = db.cursor()
    sql3 = "SELECT colname1,colname2 FROM boxplot_colname WHERE dataset = " + "'" + str(d) + "'" + " AND area = " + "'" + str(c) + "'" + " AND drug = " + "'" +str(b) + "'"
    cursor3.execute(sql3)
    data3 = cursor3.fetchall()
    colname1 = data3[0][0]
    colname2 = data3[0][1]
    db.close()
    # 获取差异分析方法
    db = pymysql.connect(**config)
    cursor9 = db.cursor()
    sql9 = "SELECT type FROM dataset_type WHERE dataset = " + "'" + str(d) + "'"
    cursor9.execute(sql9)
    data9 = cursor9.fetchall()
    db.close()
    mothod = data9[0][0]

    # 获取同源基因
    if a=="mouse":
        db = pymysql.connect(**config)
        cursor4 = db.cursor()
        sql4 = "SELECT human,mouse_ID,human_ID FROM homologene WHERE mouse = " + "'" + str(gene) + "'"
        cursor4.execute(sql4)
        data4 = cursor4.fetchall()
        if data4!=():
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
    sql5 = "SELECT `Inference Score`,`Reference Count` FROM ctd_cocaine WHERE `Gene Symbol` = " + "'" + str(homologene_human) + "' LIMIT 1"
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
    sql6 = "SELECT `Inference Score`,`Reference Count` FROM ctd_heroin WHERE `Gene Symbol` = " + "'" + str(homologene_human) + "' LIMIT 1"
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
    sql7 = "SELECT `Inference Score`,`Reference Count` FROM ctd_morphine WHERE `Gene Symbol` = " + "'" + str(homologene_human) + "' LIMIT 1"
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
    sql8 = "SELECT Organism,Interaction,InteractionActions,PubMedIDs FROM ctd_drug_gene_link WHERE `GeneSymbol` = " + "'" + str(
        homologene_human) + "'" + "AND `ChemicalName` = " + "'" + str(b) + "'"
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

    return render(request, 'search_gene.html', locals())

def geneSearch(request):
    gene1 = request.GET.get('gene')
    if gene1 == None:
        return render(request, 'geneSearch.html')
    # limma
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    sql1 = "SELECT table_name FROM dataset_type WHERE type = 'limma'"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    db.close()
    # DESeq2
    db = pymysql.connect(**config)
    cursor2 = db.cursor()
    sql2 = "SELECT table_name FROM dataset_type WHERE type = 'DESeq2'"
    cursor2.execute(sql2)
    data2 = cursor2.fetchall()
    db.close()
    table_list = []
    table_list2 = []
    for z1 in data1:
        for i1 in z1:
            db = pymysql.connect(**config)
            cursor3 = db.cursor()
            sql3 = "SELECT `Gene`,`logFC`,`P.Value`,`adj.P.Val`,`Species`,`Drug`,`Area`,`Dataset` FROM " + i1 + " WHERE Gene = " + "'" + str(gene1) + "'" + "AND `P.Value` < 0.05 And ABS(`logFC`) > log2(1.2) "
            cursor3.execute(sql3)
            data3 = cursor3.fetchall()
            db.close()
            if data3 != ():
                for i3 in data3:
                    for j2 in i3:
                        table_list2.append(j2)
                    table_list2.append(i1)
                    table_list.append(table_list2)
                    table_list2 = []

    for z2 in data2:
        for i2 in z2:
            db = pymysql.connect(**config)
            cursor4 = db.cursor()
            sql4 = "SELECT `Gene`,`log2FoldChange`,`pvalue`,`padj`,`Species`,`Drug`,`Area`,`Dataset` FROM " + i2 + " WHERE Gene = " + "'" + str(gene1) + "'" + "AND `pvalue` < 0.05 And ABS(`log2FoldChange`) > log2(1.2) "
            cursor4.execute(sql4)
            data4 = cursor4.fetchall()
            db.close()
            if data4 != ():
                for i in data4:
                    for j in i:
                        table_list2.append(j)
                    table_list2.append(i2)
                    table_list.append(table_list2)
                    table_list2 = []

    db = pymysql.connect(**config)
    cursor5 = db.cursor()
    sql5 = "SELECT table_name from information_schema.columns where table_name like 'methylation%' group by table_name;"
    cursor5.execute(sql5)
    data5 = cursor5.fetchall()
    db.close()
    methylation = []
    for i in data5:
        db = pymysql.connect(**config)
        cursor6 = db.cursor()
        sql6 = "SELECT `Probe ID`,`logFC`,`P.Value`,`deltaBeta`,`gene`,`feature`,`Drug`,`Area`,`Dataset` from " + str(i[0]) + " where gene= " + "'" + str(gene1) + "'"
        cursor6.execute(sql6)
        data6 = cursor6.fetchall()
        db.close()
        if data6 != ():
            methylation_data = 1
            for j in data6:
                methylation.append(j)
    # 生成分页器对象
    paginator = Paginator(table_list, 10)
    # 获取浏览器端请求的页码,需要设置默认值
    current_page_num = int(request.GET.get('page', 1))
    # 当前页的所有书对象
    current_page = paginator.page(current_page_num)
    # 页码列表，可迭代
    # 分页过多，需要用条件判断显示的页码
    if paginator.num_pages > 11:  # 11就是显示的固定个数
        if current_page_num - 5 < 1:  # 接近最小页码时，固定页码，否则会出现负数
            page_range = range(1, 12)
        elif current_page_num + 5 > paginator.num_pages:  # 接近最大页码时，根据最大页码限制页码数，否则会出现不存在的页码
            page_range = range(paginator.num_pages - 10, paginator.num_pages + 1)
        else:
            page_range = range(current_page_num - 5, current_page_num + 6)  # 其他情况，显示当前的挨着的几个
    else:
        page_range = paginator.page_range  # 页码总数不足时，显示全部，即不会超宽
    return render(request, 'geneSearch.html', locals())


def methylation(request):
    search_data1 = request.GET.get('a')
    search_data2 = request.GET.get('b')
    search_data3 = request.GET.get('c')
    p = request.GET.get('p')
    deltaBeta = request.GET.get('deltaBeta')

    if search_data1 == None:
        return render(request, 'methylationSearch.html')


    sql11 = "SELECT `Probe ID`,`logFC`,`AveExpr`,`t`,`P.Value`,`B`,`deltaBeta`,`gene`,`feature`,`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(
        search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + " WHERE `P.Value` < 0.05 And ABS(`deltaBeta`) > 0.05 "
    sql12 = "SELECT `Probe ID`,`logFC`,`AveExpr`,`t`,`P.Value`,`B`,`deltaBeta`,`gene`,`feature`,`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(
        search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + " WHERE `P.Value` < 0.05 And ABS(`deltaBeta`) > 0.1 "
    sql13 = "SELECT `Probe ID`,`logFC`,`AveExpr`,`t`,`P.Value`,`B`,`deltaBeta`,`gene`,`feature`,`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(
        search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + " WHERE `P.Value` < 0.05 And ABS(`deltaBeta`) > 0.2 "
    sql21 = "SELECT `Probe ID`,`logFC`,`AveExpr`,`t`,`P.Value`,`B`,`deltaBeta`,`gene`,`feature`,`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(
        search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + " WHERE `P.Value` < 0.01 And ABS(`deltaBeta`) > 0.05 "
    sql22 = "SELECT `Probe ID`,`logFC`,`AveExpr`,`t`,`P.Value`,`B`,`deltaBeta`,`gene`,`feature`,`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(
        search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + " WHERE `P.Value` < 0.01 And ABS(`deltaBeta`) > 0.1 "
    sql23 = "SELECT `Probe ID`,`logFC`,`AveExpr`,`t`,`P.Value`,`B`,`deltaBeta`,`gene`,`feature`,`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(
        search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + " WHERE `P.Value` < 0.01 And ABS(`deltaBeta`) > 0.2 "

    db = pymysql.connect(**config)
    cursor11 = db.cursor()
    cursor11.execute(sql11)
    data11 = cursor11.fetchall()
    db.close()
    n11 = len(data11)

    db = pymysql.connect(**config)
    cursor12 = db.cursor()
    cursor12.execute(sql12)
    data12 = cursor12.fetchall()
    db.close()
    n12 = len(data12)

    db = pymysql.connect(**config)
    cursor13 = db.cursor()
    cursor13.execute(sql13)
    data13 = cursor13.fetchall()
    db.close()
    n13 = len(data13)

    db = pymysql.connect(**config)
    cursor21 = db.cursor()
    cursor21.execute(sql21)
    data21 = cursor21.fetchall()
    db.close()
    n21 = len(data21)

    db = pymysql.connect(**config)
    cursor22 = db.cursor()
    cursor22.execute(sql22)
    data22 = cursor22.fetchall()
    db.close()
    n22 = len(data22)

    db = pymysql.connect(**config)
    cursor23 = db.cursor()
    cursor23.execute(sql23)
    data23 = cursor23.fetchall()
    db.close()
    n23 = len(data23)

    # ProbeLasso算法数量
    db = pymysql.connect(**config)
    cursor3 = db.cursor()
    sql3 = "SELECT * FROM m_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
        search_data3) + "_dmr_p"
    cursor3.execute(sql3)
    data3 = cursor3.fetchall()
    db.close()
    ProbeLasso_lenth = len(data3)

    # 火山图数据
    db = pymysql.connect(**config)
    cursor2 = db.cursor()
    if p:
        sql2 = "SELECT `volcano plot` FROM plot_methylation WHERE `dataset name`=" + "'" + str(
            search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_" + str(p) + "_" + str(deltaBeta)+"'"
    else:
        sql2 = "SELECT `volcano plot` FROM plot_methylation WHERE `dataset name`=" + "'" + str(
            search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_0.05_0.05'"
    cursor2.execute(sql2)
    data2 = cursor2.fetchall()
    db.close()
    if data2 != None:
        volcano_plot_base64 = data2[0][0]

    if p:
        if p == '0.05' and deltaBeta == '0.05':
            go_table_name = "m_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + "_go"
            kegg_table_name = "m_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + "_kegg"
        else:
            go_table_name = "m_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + "_go_" + str(p) + "_" + str(deltaBeta)
            kegg_table_name = "m_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + "_kegg_" + str(p) + "_" + str(deltaBeta)
    else:
        go_table_name = "m_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_go"
        kegg_table_name = "m_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
            search_data3) + "_kegg"

        # go-BP富集分析
    db = pymysql.connect(**config)
    cursor21 = db.cursor()
    sql21 = "SELECT `Size`,`Description`,`pvalue`,`ID` FROM `" + str(go_table_name) + "` WHERE ONTOLOGY = 'BP' LIMIT 10"
    cursor21.execute(sql21)
    go_bp = cursor21.fetchall()
    db.close()
    BP_list = []
    BP_list1 = []
    GO_name = []
    GO_count = []
    GO_padjust = []
    GO_list = []
    for z in go_bp:
        GO_count.append(z[0])
        GO_name.append(z[1])
        GO_padjust.append(z[2])
        for y in z:
            BP_list1.append(y)
        BP_list.append(BP_list1)
        GO_list.append(BP_list1)
        BP_list1 = []

    # go-CC富集分析
    db = pymysql.connect(**config)
    cursor22 = db.cursor()
    sql22 = "SELECT `Size`,`Description`,`pvalue`,`ID` FROM `" + str(go_table_name) + "` WHERE ONTOLOGY = 'CC' LIMIT 10"
    cursor22.execute(sql22)
    go_cc = cursor22.fetchall()
    db.close()
    CC_list = []
    CC_list1 = []
    for z in go_cc:
        GO_count.append(z[0])
        GO_name.append(z[1])
        GO_padjust.append(z[2])
        for y in z:
            CC_list1.append(y)
        CC_list.append(CC_list1)
        GO_list.append(CC_list1)
        CC_list1 = []

    # go-MF富集分析
    db = pymysql.connect(**config)
    cursor23 = db.cursor()
    sql23 = "SELECT `Size`,`Description`,`pvalue`,`ID` FROM `" + str(go_table_name) + "` WHERE ONTOLOGY = 'MF' LIMIT 10"
    cursor23.execute(sql23)
    go_mf = cursor23.fetchall()
    db.close()
    MF_list = []
    MF_list1 = []
    for z in go_mf:
        GO_count.append(z[0])
        GO_name.append(z[1])
        GO_padjust.append(z[2])
        for y in z:
            MF_list1.append(y)
        MF_list.append(MF_list1)
        GO_list.append(MF_list1)
        MF_list1 = []

    if GO_count:
        GO_count_max = max(GO_count)
        GO_count_min = min(GO_count)
    if GO_padjust:
        GO_padjust_max = max(GO_padjust)
        GO_padjust_min = min(GO_padjust)
    GO_name.reverse()

    # kegg富集分析
    db = pymysql.connect(**config)
    cursor20 = db.cursor()
    sql20 = "SELECT `Size`,`Description`,`pvalue`,`ID` FROM `" + str(kegg_table_name) + "` WHERE `pvalue` < 0.05 LIMIT 30"
    cursor20.execute(sql20)
    kegg = cursor20.fetchall()
    db.close()
    kegg_list = []
    kegg_list1 = []
    kegg_name = []
    kegg_count = []
    kegg_pvalue = []
    for z in kegg:
        kegg_count.append(z[0])
        kegg_name.append(z[1])
        kegg_pvalue.append(z[2])
        for y in z:
            kegg_list1.append(y)
        kegg_list.append(kegg_list1)
        kegg_list1 = []
    if kegg != ():
        kegg_count_max = max(kegg_count)
        kegg_count_min = min(kegg_count)
        kegg_pvalue_max = max(kegg_pvalue)
        kegg_pvalue_min = min(kegg_pvalue)
        kegg_name.reverse()

    return render(request, 'methylationSearch.html', locals())


def dmp(request):
    """ 搜索 """
    search_data1 = request.GET.get('a')
    search_data2 = request.GET.get('b')
    search_data3 = request.GET.get('c')
    p = request.GET.get('p')
    deltaBeta = request.GET.get('deltaBeta')

    search_data0 = request.GET.get('filter')

    sort = request.GET.get('sort')
    dir = request.GET.get('dir')

    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if sort == 'Pvalue':
        sotr = 'P.Value'
    elif sort == 'logFC':
        sotr = 'logFC'
    else:
        sotr = 'deltaBeta'
    if search_data0 != '' and search_data0 != None:
        if dir == 'desc':
            sql = "SELECT `Probe ID`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                  "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                  "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                  "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                  "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B`," \
                  "CASE WHEN ABS(`deltaBeta`) < 1E-8 THEN CAST(`deltaBeta` AS CHAR) ELSE CAST(`deltaBeta` AS DECIMAL(18, 8)) END AS `deltaBeta`,`gene`,`feature`," \
                  "`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + " WHERE `P.Value` < " + p + " And ABS(`deltaBeta`) > " + deltaBeta + " And (`Gene` LIKE " + "'%" + str(
                search_data0) + "%' OR `Probe ID` LIKE " + "'%" + str(
                search_data0) + "%' OR `feature` LIKE " + "'%" + str(
                search_data0) + "%' OR `cgi` LIKE " + "'%" + str(search_data0) + "%') ORDER BY CONVERT(`" + str(
                sotr) + "`, DOUBLE) DESC;"
        elif dir == 'asc':
            sql = "SELECT `Probe ID`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                  "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                  "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                  "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                  "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B`," \
                  "CASE WHEN ABS(`deltaBeta`) < 1E-8 THEN CAST(`deltaBeta` AS CHAR) ELSE CAST(`deltaBeta` AS DECIMAL(18, 8)) END AS `deltaBeta`,`gene`,`feature`," \
                  "`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + " WHERE `P.Value` < " + p + " And ABS(`deltaBeta`) > " + deltaBeta + " And (`Gene` LIKE " + "'%" + str(
                search_data0) + "%' OR `Probe ID` LIKE " + "'%" + str(
                search_data0) + "%' OR `feature` LIKE " + "'%" + str(
                search_data0) + "%' OR `cgi` LIKE " + "'%" + str(search_data0) + "%') ORDER BY CONVERT(`" + str(
                sotr) + "`, DOUBLE) ASC;"
        else:
            sql = "SELECT `Probe ID`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                  "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                  "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                  "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                  "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B`," \
                  "CASE WHEN ABS(`deltaBeta`) < 1E-8 THEN CAST(`deltaBeta` AS CHAR) ELSE CAST(`deltaBeta` AS DECIMAL(18, 8)) END AS `deltaBeta`,`gene`,`feature`," \
                  "`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + " WHERE `P.Value` < " + p + " And ABS(`deltaBeta`) > " + deltaBeta + " And (`Gene` LIKE " + "'%" + str(
                search_data0) + "%' OR `Probe ID` LIKE " + "'%" + str(
                search_data0) + "%' OR `feature` LIKE " + "'%" + str(
                search_data0) + "%' OR `cgi` LIKE " + "'%" + str(search_data0) + "%')"
    else:
        if dir == 'desc':
            sql = "SELECT `Probe ID`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                  "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                  "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                  "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                  "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B`," \
                  "CASE WHEN ABS(`deltaBeta`) < 1E-8 THEN CAST(`deltaBeta` AS CHAR) ELSE CAST(`deltaBeta` AS DECIMAL(18, 8)) END AS `deltaBeta`,`gene`,`feature`," \
                  "`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + " WHERE `P.Value` < " + p + " And ABS(`deltaBeta`) > " + deltaBeta + "ORDER BY CONVERT(`" + str(
                sotr) + "`, DOUBLE) DESC;"
        elif dir == 'asc':
            sql = "SELECT `Probe ID`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                  "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                  "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                  "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                  "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B`," \
                  "CASE WHEN ABS(`deltaBeta`) < 1E-8 THEN CAST(`deltaBeta` AS CHAR) ELSE CAST(`deltaBeta` AS DECIMAL(18, 8)) END AS `deltaBeta`,`gene`,`feature`," \
                  "`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + " WHERE `P.Value` < " + p + " And ABS(`deltaBeta`) > " + deltaBeta + "ORDER BY CONVERT(`" + str(
                sotr) + "`, DOUBLE) ASC;"
        else:
            sql = "SELECT `Probe ID`,CASE WHEN ABS(`logFC`) < 1E-8 THEN CAST(`logFC` AS CHAR) ELSE CAST(`logFC` AS DECIMAL(18, 8)) END AS `logFC`," \
                  "CASE WHEN ABS(`AveExpr`) < 1E-8 THEN CAST(`AveExpr` AS CHAR) ELSE CAST(`AveExpr` AS DECIMAL(18, 8)) END AS `AveExpr`," \
                  "CASE WHEN ABS(`t`) < 1E-8 THEN CAST(`t` AS CHAR) ELSE CAST(`t` AS DECIMAL(18, 8)) END AS `t`," \
                  "CASE WHEN ABS(`P.Value`) < 1E-8 THEN CAST(`P.Value` AS CHAR) ELSE CAST(`P.Value` AS DECIMAL(18, 8)) END AS `P.Value`," \
                  "CASE WHEN ABS(`B`) < 1E-8 THEN CAST(`B` AS CHAR) ELSE CAST(`B` AS DECIMAL(18, 8)) END AS `B`," \
                  "CASE WHEN ABS(`deltaBeta`) < 1E-8 THEN CAST(`deltaBeta` AS CHAR) ELSE CAST(`deltaBeta` AS DECIMAL(18, 8)) END AS `deltaBeta`,`gene`,`feature`," \
                  "`cgi`,`UCSC_CpG_Islands_Name` FROM methylation_" + str(
                search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + " WHERE `P.Value` < " + p + " And ABS(`deltaBeta`) > " + deltaBeta
    cursor1.execute(sql)
    data = cursor1.fetchall()
    db.close()

    # 生成分页器对象
    paginator = Paginator(data, 10)
    # 获取浏览器端请求的页码,需要设置默认值
    current_page_num = int(request.GET.get('page', 1))
    # 当前页的所有书对象
    current_page = paginator.page(current_page_num)
    # 页码列表，可迭代
    # 分页过多，需要用条件判断显示的页码
    if paginator.num_pages > 11:  # 11就是显示的固定个数
        if current_page_num - 5 < 1:  # 接近最小页码时，固定页码，否则会出现负数
            page_range = range(1, 12)
        elif current_page_num + 5 > paginator.num_pages:  # 接近最大页码时，根据最大页码限制页码数，否则会出现不存在的页码
            page_range = range(paginator.num_pages - 10, paginator.num_pages + 1)
        else:
            page_range = range(current_page_num - 5, current_page_num + 6)  # 其他情况，显示当前的挨着的几个
    else:
        page_range = paginator.page_range  # 页码总数不足时，显示全部，即不会超宽
    return render(request, 'methylation_DMP.html', locals())

def dmpgene(request):
    gene = request.GET.get('gene')
    a = request.GET.get('a')
    b = request.GET.get('b')
    c = request.GET.get('c')
    # 将小写变成大写
    c_upper = c.upper()
    # 获取同源基因
    db = pymysql.connect(**config)
    cursor = db.cursor()
    sql = "SELECT mouse,mouse_ID,human_ID FROM homologene WHERE human = " + "'" + str(gene) + "'"
    cursor.execute(sql)
    data = cursor.fetchall()
    if data != ():
        homologene_mouse = data[0][0]
        mouse_ID = data[0][1]
        human_ID = data[0][2]
    else:
        homologene_mouse = gene
    homologene_human = gene
    db.close()

    # 获取ctd-cocaine
    cocaine_list = []
    db = pymysql.connect(**config)
    cursor5 = db.cursor()
    sql5 = "SELECT `Inference Score`,`Reference Count` FROM ctd_cocaine WHERE `Gene Symbol` = " + "'" + str(homologene_human) + "' LIMIT 1"
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
    sql6 = "SELECT `Inference Score`,`Reference Count` FROM ctd_heroin WHERE `Gene Symbol` = " + "'" + str(homologene_human) + "' LIMIT 1"
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
    sql7 = "SELECT `Inference Score`,`Reference Count` FROM ctd_morphine WHERE `Gene Symbol` = " + "'" + str(homologene_human) + "' LIMIT 1"
    cursor7.execute(sql7)
    data7 = cursor7.fetchall()
    for i in data7:
        morphine_list.append(float(i[0]))
        morphine_list.append(float(i[1]))
    db.close()

    # 获取ctd-amphetamine
    amphetamine_list = []
    db = pymysql.connect(**config)
    cursor4 = db.cursor()
    sql4 = "SELECT `Inference Score`,`Reference Count` FROM ctd_amphetamine WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor4.execute(sql4)
    data4 = cursor4.fetchall()
    for i in data4:
        amphetamine_list.append(float(i[0]))
        amphetamine_list.append(float(i[1]))
    db.close()

    # 获取ctd-cannabis
    cannabis_list = []
    db = pymysql.connect(**config)
    cursor3 = db.cursor()
    sql3 = "SELECT `Inference Score`,`Reference Count` FROM ctd_cannabis WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor3.execute(sql3)
    data3 = cursor3.fetchall()
    for i in data3:
        cannabis_list.append(float(i[0]))
        cannabis_list.append(float(i[1]))
    db.close()

    # 获取ctd-ethanol
    ethanol_list = []
    db = pymysql.connect(**config)
    cursor2 = db.cursor()
    sql2 = "SELECT `Inference Score`,`Reference Count` FROM ctd_ethanol WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor2.execute(sql2)
    data2 = cursor2.fetchall()
    for i in data2:
        ethanol_list.append(float(i[0]))
        ethanol_list.append(float(i[1]))
    db.close()

    # 获取ctd-nicotine
    nicotine_list = []
    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    sql1 = "SELECT `Inference Score`,`Reference Count` FROM ctd_nicotine WHERE `Gene Symbol` = " + "'" + str(
        homologene_human) + "' LIMIT 1"
    cursor1.execute(sql1)
    data1 = cursor1.fetchall()
    for i in data1:
        nicotine_list.append(float(i[0]))
        nicotine_list.append(float(i[1]))
    db.close()

    # 获取ctd中药物与基因的文献信息
    db = pymysql.connect(**config)
    cursor8 = db.cursor()
    sql8 = "SELECT Organism,Interaction,InteractionActions,PubMedIDs FROM ctd_drug_gene_link WHERE `GeneSymbol` = " + "'" + str(
        homologene_human) + "'" + "AND `ChemicalName` = " + "'" + str(a) + "'"
    cursor8.execute(sql8)
    data8 = cursor8.fetchall()
    db.close()

    return render(request, 'methylation_DMP_gene.html', locals())

def dmr(request):
    """ 搜索 """
    search_data1 = request.GET.get('a')
    search_data2 = request.GET.get('b')
    search_data3 = request.GET.get('c')
    method = request.GET.get('method')

    search_data0 = request.GET.get('filter')

    db = pymysql.connect(**config)
    cursor1 = db.cursor()
    if search_data0 != '' and search_data0 != None:
        if method == 'Bumhunter':
            sql1 = "SELECT * FROM m_" + str(
                search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_dmr_b WHERE `BumphunterDMR.seqnames` LIKE " + "'%" + str(
                search_data0) + "%'"
        else:
            sql1 = "SELECT `ProbeLassoDMR.dmrChrom`,`ProbeLassoDMR.dmrP`,`ProbeLassoDMR.dmrpRank`,`ProbeLassoDMR.dmrStart`,`ProbeLassoDMR.dmrEnd`,`ProbeLassoDMR.dmrSize`,`ProbeLassoDMR.dmrCoreStart`,`ProbeLassoDMR.dmrCoreEnd`,`ProbeLassoDMR.dmrCoreSize`,`ProbeLassoDMR.ensemblID`,`ProbeLassoDMR.geneSymbol` FROM m_" + str(
                search_data1) + "_" + str(search_data2) + "_" + str(
                search_data3) + "_dmr_p WHERE `ProbeLassoDMR.seqnames` LIKE " + "'%" + str(
                search_data0) + "%' OR`ProbeLassoDMR.geneSymbol` LIKE " + "'%" + str(
                search_data0) + "%'"
    else:
        if method == 'Bumhunter':
            sql1 = "SELECT * FROM m_" + str(
                search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_dmr_b"
        else:
            sql1 = "SELECT `ProbeLassoDMR.dmrChrom`,`ProbeLassoDMR.dmrP`,`ProbeLassoDMR.dmrpRank`,`ProbeLassoDMR.dmrStart`,`ProbeLassoDMR.dmrEnd`,`ProbeLassoDMR.dmrSize`,`ProbeLassoDMR.dmrCoreStart`,`ProbeLassoDMR.dmrCoreEnd`,`ProbeLassoDMR.dmrCoreSize`,`ProbeLassoDMR.ensemblID`,`ProbeLassoDMR.geneSymbol` FROM m_" + str(
                search_data1) + "_" + str(search_data2) + "_" + str(search_data3) + "_dmr_p"

    cursor1.execute(sql1)
    data = cursor1.fetchall()
    db.close()

    # 生成分页器对象
    paginator = Paginator(data, 10)
    # 获取浏览器端请求的页码,需要设置默认值
    current_page_num = int(request.GET.get('page', 1))
    # 当前页的所有书对象
    current_page = paginator.page(current_page_num)
    # 页码列表，可迭代
    # 分页过多，需要用条件判断显示的页码
    if paginator.num_pages > 11:  # 11就是显示的固定个数
        if current_page_num - 5 < 1:  # 接近最小页码时，固定页码，否则会出现负数
            page_range = range(1, 12)
        elif current_page_num + 5 > paginator.num_pages:  # 接近最大页码时，根据最大页码限制页码数，否则会出现不存在的页码
            page_range = range(paginator.num_pages - 10, paginator.num_pages + 1)
        else:
            page_range = range(current_page_num - 5, current_page_num + 6)  # 其他情况，显示当前的挨着的几个
    else:
        page_range = paginator.page_range  # 页码总数不足时，显示全部，即不会超宽
    return render(request, 'methylation_DMR.html', locals())