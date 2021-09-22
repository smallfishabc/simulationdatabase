import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib
import numpy as np

# matplotlib.rcParams['svg.fonttype'] = 'none'
# proteinlist=['E1A','Ash1','p53','puma_wild']
# data = pd.read_csv('tryyourbest.csv')
def plotScatter(data,name,metric,proteins='.*',colorBy='idx',legend=0,xlabel='',ylabel='',ylimit=''):
    fix,ax=plt.subplots(1,len(metric),figsize=[3*len(metric),5])
    plt.tight_layout()
    B=data[data['Protein name'].str.match(proteins)].reset_index()
    B['idx']=B['Protein name']
    print(B)
    plt.figure(figsize=(20,20))
    order=1
    colortitle=['Attractive solution','Repulsive solution']
    colorref=['b','r']
    standard=['p53','E1A','puma_wildfull','Ash1']
    Highlight=data[data['Protein name'].isin(standard)]
    print(Highlight)
    for j in range(len(name)):
        for i in range(len(metric)):
           # sns.scatterplot(x=B[metric[i]], y=B[name[j]], data=B, color=colorref[order])
            #plt.scatter(x=B[metric[i]], y=B[name[j]],s=B['Length'],color=colorref[order],alpha=0.25,linewidths=15,label=colortitle[order])
            plt.scatter(x=B[metric[i]]-B['Allpink(0)'], y=B[name[j]],s=B['Length'],color='black',alpha=0.25,linewidths=15,label=colortitle[order])
            sns.regplot(x=metric[i],y=name[j],data=B,color='k')
            coeff=np.corrcoef(B[metric[i]],B[name[j]])
            print(coeff)
            if order==0:
                aaa='stdChi(3-0)'
            else:
                aaa='stdChi(-3-0)'
            #plt.errorbar(x=B[metric[i]], y=B[name[j]],yerr=B[aaa],color=colorref[order],alpha=0.15,fmt='none',label=None)
            plt.errorbar(x=B[metric[i]]-B['Allpink(0)'], y=B[name[j]],yerr=B['stdChi(0)'],color='black',alpha=0.15,fmt='none',label=None)
#            plt.scatter(x=Highlight[metric[i]], y=Highlight[name[j]],marker='s',linewidths=20)
            ax=plt.gca()


def plot_modification(fig,ax,xlabel,ylabel,ylimit):
    ax.tick_params(direction='out', length=12, width=2)
    plt.xticks([-0.5,0,0.5,1],fontsize=40)
    plt.yticks(fontsize=40)
    plt.xticks(fontsize=40)
    plt.xlabel(xlabel, fontsize=40)
    plt.ylabel(ylabel, fontsize=60)
    plt.ylim(-1.5,ylimit)
    plt.legend(fontsize=40)

def custom_plot():
    return