"""
URL configuration for database project.

The `urlpatterns` list routes URLs to views. For more information please see:
    https://docs.djangoproject.com/en/4.2/topics/http/urls/
Examples:
Function views
    1. Add an import:  from my_app import views
    2. Add a URL to urlpatterns:  path('', views.home, name='home')
Class-based views
    1. Add an import:  from other_app.views import Home
    2. Add a URL to urlpatterns:  path('', Home.as_view(), name='home')
Including another URLconf
    1. Import the include() function: from django.urls import include, path
    2. Add a URL to urlpatterns:  path('blog/', include('blog.urls'))
"""
from django.urls import path, re_path
from drug_addiction.views import search, index, views

urlpatterns = [
    path('', index.index),
    path('datasetSearch/', search.datasetSearch),
    # path('centrality_select', search.centrality_select),
    path('index/', index.index),
    path('search/', search.AreaView.as_view()),
    path('getInfo/', search.getInfo),
    path('datasetSearch/DEGs/gene/', search.gene),
    path('geneSearch/', search.geneSearch),
    path('datasetSearch/DEGs/', search.degs),
    path('help/', views.help),
    path('methylationSearch/', search.methylation),
    path('getMethylationInfo/', search.getMethylationInfo),
    path('methylationSearch/DMP/', search.dmp),
    path('methylationSearch/DMP/gene/', search.dmpgene),
    path('methylationSearch/DMR/', search.dmr),
    path('methylationSearch/DMR/gene/', search.dmpgene),
    path('download/', views.download),
    path('methylationdownload/', views.methylationdownload),
    path('download/degs', views.download1, name='download1'),
    path('download/ppi', views.download2, name='download2'),
    path('download/centrality', views.download3, name='download3'),
    path('download/go', views.download4, name='download4'),
    path('download/kegg', views.download5, name='download5'),
    path('methylationdownload/dmg', views.download6, name='download6'),
    path('methylationdownload/dmr/Bumphunter', views.download7, name='download7'),
    path('methylationdownload/dmr/ProbeLasso', views.download8, name='download8'),
    path('methylationdownload/go', views.download9, name='download9'),
    path('methylationdownload/kegg', views.download10, name='download10'),
    path('contact/', views.contact),

    path('upload/', index.upload),
    path('analyseResult/', index.success_view, name='success'),
    path('analyseResult/allResult', index.download, name='download'),
    path('analyseResult/keyGene/', index.keygene),
]
