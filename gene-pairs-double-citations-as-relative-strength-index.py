import os
import glob
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from math import log
import pickle
from selenium import webdriver
from selenium.webdriver.common.keys import Keys
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.common.exceptions import NoSuchElementException

def createOrganSpecificGenesFile(organ,filename):
    "scrapes a list of genes from TheHUmanProteeinAtlas according to specified organ and saves as .txt file"
    
    driver = webdriver.Chrome(ChromeDriverManager().install()) 
    driver.get("https://www.proteinatlas.org/search/tissue_category_rna%3A"+organ+"%3BTissue+enriched%2CGroup+enriched%2CTissue+enhanced+AND+sort_by%3Atissue+specific+score/2")
    genes=[]
    for page in range(10):
        raw_genes = driver.find_elements_by_class_name("tda");
        for i in range(len(raw_genes)):
            if i%3 ==0:
                genes.append(raw_genes[i].text);
        next_page = driver.find_element_by_partial_link_text('next')
        next_page.click()
    print("found "+len(genes)+" genes")
    textfile=open(organ+'_'+filename+'.txt','w')
    for gene in genes:
        textfile.write(gene + '\n')
    textfile.close()
    return genes

def createSingleCitation(genes,organ):
    "finds the number of abstracts in pubmed that specify the gene in their abstract"
    "returns a df of genes and their respective number of abstract citations and saves .csv"
    
    # driver = webdriver.Chrome(ChromeDriverManager().install()) 
    single_citations = pd.DataFrame(index=genelist, columns = ['Genes','citations #'])
    driver = webdriver.Chrome(ChromeDriverManager().install())
    for i,gene in enumerate(genelist):
        driver.get('https://pubmed.ncbi.nlm.nih.gov/?term='+gene+'%5BText+Word%5D&sort=')
    # driver.get("https://pubmed.ncbi.nlm.nih.gov/?term="+gene)
        try:
            publication = driver.find_element(By.ID, 'search-results').text.split(' ')[0]
            if publication=='No':
                single_citations['citations #'][gene] = 0
            else:
                single_citations['citations #'][gene] = int(publication.replace(',',''))
        except NoSuchElementException:
            single_citations['citations #'][gene] = 1
        if i%25 ==0:
            print(single_citations[i-25:i])
    single_citation_sorted = single_citations.sort_values('citations #', axis = 0,ascending = False).reset_index()
    single_citation_sorted.to_csv(organ+'_Single_Citations.csv')
    return(single_citation_sorted)

def createDoubleCitations(genes,single_citation_sorted,organ):
    "finds the number of abstracts that have two specified genes in pubmed, for each two genes in a list"
    "creates a mtrix of double citations are saves a csv file"
    
    driver = webdriver.Chrome(ChromeDriverManager().install()) 
    n_genes = len(genes)
    cocit_mat = np.zeros((n_genes,n_genes))-1
    for i in range(n_genes):
        geneA = single_citation_sorted['Genes'][i]
        for j in range(i):
            geneB = single_citation_sorted['Genes'][j]
            driver.get('https://pubmed.ncbi.nlm.nih.gov/?term='+geneA+'+'+geneB+'%5BText+Word%5D&sort=')
            try:
                publication = driver.find_element(By.ID, 'search-results').text.split(' ')[0]
                if publication == 'No':
                    cocit_mat[i,j] = 0
                else:
                    cocit_mat[i,j] = int(publication.replace(',', ''))
            except NoSuchElementException:
                cocit_mat[i,j] = 1
        print(i,' ',geneA,' max co-pubs: ',np.max(cocit_mat[i,:]))
        # if i%25 ==0:
    print(cocit_mat[:i, :i])
    print('***')
    cocit_mat.tofile('Cocitation_Matrix'+str(gene_number_to_start+batch)+'.csv',sep=',',format='%10.1f')
    return cocit_mat

def calcDoubleCitationStrength(cocit_mat,single_citation_sorted,total_number_of_abstracts):
    "calculates the double citation strength. i.e. the nubmer of co-publication the gene pair has in relation to
    "the value of its single citations aka normalised_strength"
    single_citation_np= = single_citation_sorted['citations #'].to_numpy()
    np.fill_diagonal(cocit_mat,single_citation_np)
    cocit_strength = np.zeros_like(cocit_mat)
    cocit_stength_distribution =[]
    for i in range(len(cocit_strength)):
        for j in range(i,len(cocit_strength)):
            if i!=j:
                cocit_strength[i,j] = cocit_mat[i,j]*total_number_of_abstracts/(cocit_mat[i,i]*cocit_mat[j,j])
                if cocit_strength[i,j]>0:
                    cocit_stength_distribution.append(cocit_strength[i,j])
                cocit_strength[j,i] = cocit_strength[i,j]
    np.save('cocitation_matrix_normalised',cocit_strength,allow_pickle=True, fix_imports=True)
    return cocit_strength, cocit_stength_distribution
   
def heatmapLogValues(cocit_strength):
    ax = sns.heatmap(np.log(cocit_strength+0.000001), cmap="magma", vmin = -1, vmax= 6, square=True, linewidth=0.00005)
    plt.show()
    return
    
def binIndicesBasedOnStrength(cocit_stength_distribution):
    "provides the normalised_strngeth values that divide the log distribution of values into 10 equal-volume bins"
    "these values wil be used to give relative strength for each potential edge in the network, when randomly"
    "picking edges to complete build the network"
    fig, ax =plt.subplots(1,2)
    sns.set_theme(style="whitegrid")
    palette1 = sns.color_palette("tab20b")
    ax[0]=sns.histplot(cocit_stength_distribution,edgecolor=".3",color=palette1[1],linewidth=.5,log_scale=True,kde=True,ax=ax[0])
    ax[1]=sns.histplot(cocit_stength_distribution,edgecolor=".3",color=palette1[1],linewidth=.5,log_scale=True,kde=True,ax=ax[1])
    ax[0].set(xticks=[1,10,100,1000,10000,100000,1000000])
    ax[1].set(xticks=[1,10,100,1000,10000,100000,1000000])
    cocit_stength_distribution.sort()

    max_bin = np.log(np.max(cocit_stength_distribution))
    bin_indices = list(np.rint(np.exp(np.linspace(0, max_bin, num=10))))
    print(bin_indices)
    bin_indices = [int(i) for i in bin_indices]

    for i in bin_indices[:-1]:
        ax[1].axvline(i, 0, 1, c='yellow',linewidth=0.75)
    fig.tight_layout()
    plt.show()
    with open("bin_values_plink", "wb") as fp:
        pickle.dump(bin_indices, fp)
    return bin_indices
    
if __name__ == '__main__':
    #the script takes as input an organ, finds the relavent genes from the human protein atlas db that are expressed"
    #and finds the relative number each of the genes in the gene-pair in the list are published together in the same abstract"
    #these numbers will then be used as the relative chance each gene-pair will be used as an edge in a randommally contructed 
    #gene network. i.e. genes that are found together in abstracts in high frequency will have a better chance to construct the 
    # network as it grows"
    
    total_number_of_abstracts = 20000000
    organ = 'liver'
    genes = createOrganSpecificGenesFile(organ,filename)    
    single_citation_sorted = createSingleCitation(genes,organ)
    doouble_citation_matrix = createDoubleCitations(genes,single_citation_sorted,organ)
    doouble_citation_strength, strength_distribution = calcDoubleCitationStrength(cocit_mat,single_citation_sorted,total_number_of_abstracts)
    heatmapLogValues(doouble_citation_strength)
    bin_values = binIndicesBasedOnStrength(strength_distribution)
    print(bin_values)
    print("these values of the normalised strength matrix divide the distribution into 10 equally-populated bins. The bins will provide the relative chance to pick a gene pair as an edge in the network. i.e gene pairs from the fifth bin will be represented five times in the database, while pairs from the first bin will be represented only once an so on.")
    
    
    
    
