
import pandas as pd
import os
import subprocess



##Gene LOC file preparation
df=pd.read_csv("/home/jjohn41/Softwares/Resourses/AUG_MAGMA/Gene_LocFiles/NCBI37.3.gene.loc",sep="\t",header=None)
#mhc=df[ (df[1].isin(["X","Y"])) | ((df[1]=="6")  & (df[2]>=28477797) & (df[2]<=33448354))]
mhc=df[ (df[1].isin(["X","Y"])) | ((df[1]=="6")  & (df[2]>=28000000) & (df[2]<=30000000))]



###msigdb.v2023.1.Hs.symbols.gmt
os.system("cp /home/jjohn41/Softwares/Resourses/AUG_MAGMA/GenneSets/singh_et_al/singh_et_al_magma_geneset.v2.gmt singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY.gmt")
os.system("sed -i 's/,//g' singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY.gmt")

for gene in list(mhc[5]):
    os.system(f"sed -i 's/\<{gene}\>//g' singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY.gmt")

os.system(f'cat singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY.gmt| tr -s "\t" > test')
os.system(f"mv test singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY.gmt")



###################------------------------------------------------------------------------

os.system('''awk -F"\t" '{print  $1}' singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY.gmt  >msigdb.v2023.1.Hs_G2CDB_CHD8_Lisst.txt  ''')


pathdf=pd.read_csv("msigdb.v2023.1.Hs_G2CDB_CHD8_Lisst.txt",header=None)
file="singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY.gmt"

os.system(f"cut  -f1,2 {file} > {file}_first_two.tsv")
os.system(f"cut  --complement -f1,2 {file} > {file}_no_first_two.tsv")
os.system(f'cat {file}_no_first_two.tsv |tr "\t" "," > test.csv')
os.system(f'mv test.csv {file}_no_first_two.tsv ')
os.system(f'paste -d "\t" {file}_first_two.tsv {file}_no_first_two.tsv | tr " "  "\t"  > {file}.tsv')
os.system(f"rm *first_two.tsv")
df=pd.read_csv(f"{file}.tsv",sep="\t",header=None)
df=df[~df[2].isna()]
df['NumberofGenes']=df[2].apply(lambda x: len(x.split(",")))
df[[0,1,2]].to_csv(f"{file}",sep="\t",header=None,index=None)
df[[0,1,2,"NumberofGenes"]].to_csv(f"{file}.tsv",sep="\t",index=None)
os.system(f'''sed -i 's/,/\t/g' {file}''''')
os.system(f"rm  {file}.tsv")
merged=pd.merge(df,pathdf,on=0)
merged[[0,1,2]].to_csv(f"{file[:-4]}_Selected.gmt",sep="\t",header=None,index=None)
merged[[0,1,2,'NumberofGenes']].to_csv(f"{file[:-4]}_Selected.gmt.tsv",sep="\t",index=None)
os.system(f'''sed -i 's/,/\t/g' {file[:-4]}_Selected.gmt''''')
for number in [5,10,15,20]:
    test3=df[df['NumberofGenes']>=number]
    test3[[0,1,2]].to_csv(f"{file[:-4]}_Minimum{number}.gmt",sep="\t",header=None,index=None)
    test3.to_csv(f"{file[:-4]}_Minimum{number}.gmt.tsv",sep="\t",header=None,index=None)
    os.system(f'''sed -i 's/,/\t/g' {file[:-4]}_Minimum{number}.gmt''''')
    merged=pd.merge(test3,pathdf,on=0)
    merged[[0,1,2]].to_csv(f"{file[:-4]}_Minimum{number}_Selected.gmt",sep="\t",header=None,index=None)
    merged.to_csv(f"{file[:-4]}_Minimum{number}_Selected.gmt.tsv",sep="\t",header=None,index=None)
    os.system(f'''sed -i 's/,/\t/g' {file[:-4]}_Minimum{number}_Selected.gmt''''')





os.system("cp /home/jjohn41/Softwares/Resourses/AUG_MAGMA/GenneSets/singh_et_al/singh_et_al_magma_geneset.v2.gmt singh_et_al_magma_geneset.v2.symbols.gmt")
os.system("sed -i 's/,//g' singh_et_al_magma_geneset.v2.symbols.gmt")
os.system('''awk -F"\t" '{print  $1}'  singh_et_al_magma_geneset.v2.symbols.gmt  >msigdb.v2023.1.Hs_G2CDB_CHD8_Lisst.txt  ''')


pathdf=pd.read_csv("msigdb.v2023.1.Hs_G2CDB_CHD8_Lisst.txt",header=None)
file="singh_et_al_magma_geneset.v2.symbols.gmt"

os.system(f"cut  -f1,2 {file} > {file}_first_two.tsv")
os.system(f"cut  --complement -f1,2 {file} > {file}_no_first_two.tsv")
os.system(f'cat {file}_no_first_two.tsv |tr "\t" "," > test.csv')
os.system(f'mv test.csv {file}_no_first_two.tsv ')
os.system(f'paste -d "\t" {file}_first_two.tsv {file}_no_first_two.tsv | tr " "  "\t"  > {file}.tsv')
os.system(f"rm *first_two.tsv")
df=pd.read_csv(f"{file}.tsv",sep="\t",header=None)
df=df[~df[2].isna()]
df['NumberofGenes']=df[2].apply(lambda x: len(x.split(",")))
df[[0,1,2]].to_csv(f"{file}",sep="\t",header=None,index=None)
df[[0,1,2,"NumberofGenes"]].to_csv(f"{file}.tsv",sep="\t",index=None)
os.system(f'''sed -i 's/,/\t/g' {file}''''')
os.system(f"rm  {file}.tsv")
merged=pd.merge(df,pathdf,on=0)
merged[[0,1,2]].to_csv(f"{file[:-4]}_Selected.gmt",sep="\t",header=None,index=None)
merged[[0,1,2,'NumberofGenes']].to_csv(f"{file[:-4]}_Selected.gmt.tsv",sep="\t",index=None)
os.system(f'''sed -i 's/,/\t/g' {file[:-4]}_Selected.gmt''''')
for number in [5,10,15,20]:
    test3=df[df['NumberofGenes']>=number]
    test3[[0,1,2]].to_csv(f"{file[:-4]}_Minimum{number}.gmt",sep="\t",header=None,index=None)
    test3.to_csv(f"{file[:-4]}_Minimum{number}.gmt.tsv",sep="\t",header=None,index=None)
    os.system(f'''sed -i 's/,/\t/g' {file[:-4]}_Minimum{number}.gmt''''')
    merged=pd.merge(test3,pathdf,on=0)
    merged[[0,1,2]].to_csv(f"{file[:-4]}_Minimum{number}_Selected.gmt",sep="\t",header=None,index=None)
    merged.to_csv(f"{file[:-4]}_Minimum{number}_Selected.gmt.tsv",sep="\t",header=None,index=None)
    os.system(f'''sed -i 's/,/\t/g' {file[:-4]}_Minimum{number}_Selected.gmt''''')





os.sysytem('''egrep "Rbfox1|http" /home/jjohn41/Softwares/Resourses/AUG_MAGMA/GenneSets/NoCHRXY/msigdb.v2023.1.Hs_G2CDB_CHD8_AJHG.symbols_No_ChrXY_Minimum10_Selected.gmt \
                 > msigdb.v2023.1.Hs_G2CDB_CHD8_AJHG.symbols_No_ChrXY_Minimum10_Selected_toInclude.gmt''')

os.sysytem('''awk -F"\t" '{print $0}' /home/jjohn41/Softwares/Resourses/AUG_MAGMA/GenneSets/singh_et_al/singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY_Minimum10_Selected.gmt \
         |egrep -v "GOBP|GOCC|GOMF|KEGG|REACTOME|scz_gwas_genes|h.all:HALLMARK" >  singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY_Minimum10_Selected_toInclude.gmt''')

os.sysytem('''cat msigdb.v2023.1.Hs_G2CDB_CHD8_AJHG.symbols_No_ChrXY_Minimum10_Selected_toInclude.gmt  singh_et_al_magma_geneset.v2.symbols_NoMHC_ChrXY_Minimum10_Selected_toInclude.gmt >msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoChrXY_Minimum10_Selected.gmt''')




file="msigdb.v2023.1.Hs_Rbfox1_singh_custom.symbols_NoChrXY_Minimum10_Selected.gmt"

os.system(f"cut  -f1,2 {file} > {file}_first_two.tsv")
os.system(f"cut  --complement -f1,2 {file} > {file}_no_first_two.tsv")
os.system(f'cat {file}_no_first_two.tsv |tr "\t" "," > test.csv')
os.system(f'mv test.csv {file}_no_first_two.tsv ')
os.system(f'paste -d "\t" {file}_first_two.tsv {file}_no_first_two.tsv | tr " "  "\t"  > {file}.tsv')
os.system(f"rm *first_two.tsv")
df=pd.read_csv(f"{file}.tsv",sep="\t",header=None)
df=df[~df[2].isna()]
df.to_csv(f"{file[:-4]}_Minimum10.gmt.tsv",sep="\t",header=None,index=None)
