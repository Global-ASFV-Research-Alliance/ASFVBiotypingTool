#TheActualProcess
import pandas as pd
import glob
import os
from Bio import SeqIO
import argparse

argParser = argparse.ArgumentParser()
argParser.add_argument("-q", "--query", type= str, help="Input Query Fasta File Location")
argParser.add_argument("-f", "--fastabank", type= str, help="Input Fasta Folder Location", default=".\MatrixBank\\")
argParser.add_argument("-g", "--genomes", type= str, help="Input Genome CSV Location", default="./Genomes.csv")
args = argParser.parse_args()

#queryfile = 'Tests\\BioTypeRefTest\\AY261364.fa' #What they are putting in
#FastaBank =  "MatrixBank\\" #"MatrixBank\\"
#GenomeCsv = pd.read_csv('Genomes.csv', index_col='Accession')

def FolderPathFixer(FolderPath):
    """
    An internal function used to fix folder paths.
    
    Parameters
    ---
    A path to a folder.
    
    Returns
    ---
    The folder path, if not made, will be made. Additionally, a "\\" will be added at the end to allow more items to the end.
    """
    import os
    os.makedirs(FolderPath, exist_ok=True)
    FolderPath = FolderPath + "\\"
    return FolderPath.replace("\\\\","\\")

blastoutputbank = FolderPathFixer("./BlastOutput")
fastabank = FolderPathFixer(args.fastabank)

def BiotypingWithBlast(queryfile, blastoutputbank, FastaBank, GenomeCsv, SaveTempFiles = False):
    
    for record in SeqIO.parse(queryfile, "fasta"):
        ok = r'ACTGUWSMKRYBDHVN*'
        Type = all(c in ok for c in str(record.seq))
    import glob
    import os
    my_files = glob.glob(blastoutputbank + '*')
    for file in my_files:
        os.remove(file)

    my_files = glob.glob(FastaBank + '*')
    for file in my_files:
        filename = file.replace(FastaBank, '').replace('.fasta', '')
        outputlocation = blastoutputbank + filename + '.txt'
        #print(filename)
        import os
        if Type == True:
            from Bio.Blast.Applications import NcbiblastxCommandline
            blastcmdline = NcbiblastxCommandline(cmd='blastx', query = queryfile, subject= file, max_hsps = 1, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
            os.system(str(blastcmdline))
        if Type == False:
            from Bio.Blast.Applications import NcbiblastpCommandline
            blastcmdline = NcbiblastpCommandline(cmd='blastp', query = queryfile, subject= file, max_hsps = 1, out = outputlocation, outfmt = "6 qseqid sseqid evalue pident bitscore length qlen slen qstart qend sstart send qseq sseq", evalue=0.001)
            os.system(str(blastcmdline))

    def BlastBankToPidentMatrix(blastoutputbank):
        import glob
        import pandas as pd
        my_files = glob.glob(blastoutputbank + '*')
        DF = pd.DataFrame()
        DF2 = pd.DataFrame()
        DF3 = pd.DataFrame()
        for file in my_files:
            filename = file.replace(blastoutputbank, '').replace('.txt', '')
            prediction = pd.read_csv(file, sep = '\t', names=('qseqid', 'sseqid', 'evalue', 'pident', 'bitscore', 'length', 'qlen', 'slen', 'qstart', 'qend', 'sstart', 'send', 'qseq', 'sseq'))
            for name, pred in prediction.iterrows():
                DF.loc[pred['sseqid'],[filename]] = pred['pident']/100
                if pred['length']/pred['slen'] > 1:
                    Ratio = 1
                else:
                    Ratio = pred['length']/pred['slen']
                DF2.loc[pred['sseqid'],[filename]] = Ratio
                DF3.loc[pred['sseqid'],[filename]] = ((pred['pident']/100)**2)*Ratio
        return DF, DF2, DF3

    BM_Pident, BM_Ratio, BM_PidentRatio = BlastBankToPidentMatrix(blastoutputbank)
    
    if SaveTempFiles == True:
        BM_Pident.to_csv("TempPident.csv")
        BM_Ratio.to_csv("TempRatio.csv")
        BM_PidentRatio.to_csv("TempPidentRatio.csv")
    
    BM_Pident_Mean = BM_Pident.mean(axis=1, numeric_only = True)
    BM_Ratio_Mean = BM_Ratio.mean(axis=1, numeric_only = True)
    BM_PidentRatio_Mean = BM_PidentRatio.mean(axis=1, numeric_only = True)
    BM_Count = BM_PidentRatio.count(axis=1, numeric_only= True)

    Score = (BM_PidentRatio_Mean*(BM_Count**0.5)).sort_values(ascending=False)
    
    if SaveTempFiles == True:
        Score.to_csv("TempScore.csv")
    
    B646MatchPercentage = BM_Pident['B646L'][GenomeCsv.loc[Score.index[0]].name]*100

    if B646MatchPercentage >= 95.0:
        if B646MatchPercentage < 98.0:
            print("Warning! B646 Match Percentage is beneath 0.98. This could mean there is potential poor sequencing")
            + "\n"
        print(
            "This software finds the closest matching Biotype through a gene-by-gene match, meaning it matches the data submitted against all ASFV genes in our database."
            + "\n"
            + "Predicted Biotype: " + str(GenomeCsv.loc[Score.index[0]]['Biotype'])
            + "\n"
            + "Weighted Percent Gene-By-Gene Match (Closest in Biotype): " + str(BM_PidentRatio_Mean[Score.index[0]]*100) 
            + "\n"
            + "Percent B646L (P72) Match (Closest in Biotype): " + str(B646MatchPercentage)
            + "\n"
            + "Genes Found (Closest in Biotype): " + str(BM_Count[Score.index[0]]) 
            )
        if BM_PidentRatio_Mean[Score.index[0]] <= 0.975:
            print("Based on weighted percent identity, the genome exhibited a" + str(BM_PidentRatio_Mean[Score.index[0]]*100)% +  "similarity to its nearest match within the specified Biotype. This percent identity is below our recommended cutoff of 97.5%. Nevertheless, this outcome does not imply that the genome signifies a novel Biotype, however this genome should be analyzed using the complete algorithm described in 'Reclassification of ASFV into 7 Biotypes Using Unsupervised Machine Learning' (Dinhobl et al 2023). Please contact administration for assistance.")
    else: 
        print("Warning! B646 match percentage is " + str(B646MatchPercentage) + " which is very low. Ensure proper sequencing and cleaning before submitting.")       

BiotypingWithBlast(queryfile = args.query, blastoutputbank =  blastoutputbank, FastaBank = fastabank, GenomeCsv = pd.read_csv(args.genomes, index_col='Accession'), SaveTempFiles= True)