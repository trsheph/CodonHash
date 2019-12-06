import PIL
from PIL import Image
import imagehash
import statistics
import numpy as np
from scipy.fftpack import dct
#
# codonDict [dictionary]
# Ordered dictionary of codons by usage with value assigned equidistant between 0-255
#
codonDict={'TAG': 0, 'TGA': 4, 'AGG': 8, 'AGA': 12, 'TAA': 16, 'CGA': 20, 'CTA': 24, 'ATA': 28,
'TGT': 32, 'CCC': 36, 'CGG': 40, 'TGC': 44, 'CCT': 48, 'ACA': 52, 'TCA': 56, 'GGA': 60, 'TCT': 64,
'CCA': 68, 'TCC': 72, 'AGT': 76, 'ACT': 80, 'TCG': 84, 'CAC': 88, 'AAG': 92, 'GTA': 96, 'CTT': 100,
'CTC': 104, 'GGG': 108, 'TAC': 112, 'CAT': 116, 'TTG': 120, 'TTA': 124, 'ACG': 128, 'GTC': 132, 'GCT': 136,
'CAA': 140, 'TGG': 144, 'AGC': 148, 'TAT': 152, 'TTC': 156, 'AAT': 160, 'GAG': 164, 'GTT': 168, 'GAC': 172,
'GCA': 176, 'CGT': 180, 'AAC': 184, 'CGC': 188, 'TTT': 192, 'CCG': 196, 'ACC': 200, 'GGT': 204, 'ATC': 208,
'GCC': 212, 'GTG': 216, 'ATG': 220, 'CAG': 224, 'GGC': 228, 'ATT': 232, 'GAT': 236, 'AAA': 240, 'GCG': 244,
'GAA': 248, 'CTG': 252}
#
def hammingWrite(hashDB):
    genes=[]
    hashs=[]
    for i in hashDB:
        genes.append(i)
        hashs.append(hashDB[i])
    k=0
    kprev=33
    fOut=open('HammingClusters.txt', 'w')
    kalso=[]
    for l in range(len(hashs)):
        kprev=33
        k=0
        closest=0
        kalso=[]
        for i in range(len(hashs)):
            if i!=l:
                for j in range(32):
                    k=k+abs(int(hashs[l][j])-int(hashs[i][j]))
                if k<kprev:
                    closest=i
                    kprev=k
                    kalso=[genes[i]]
                elif k==kprev:
                    kalso.append(genes[i])
            k=0
        if kprev<6:
            fOut.write(str(genes[l])+"\t"+str(kprev)+"\t"+str(kalso)+"\n");
    fOut.close()
    return 1;

def loadGenomeData(fName):
    fIn=open(fName, 'r');
    tmpList=[];
    genes=[];
    genome=[];
    genesDB={}
    for line in fIn:
        tmpList.append(line);
    fIn.close()
    for i in range(int(len(tmpList)/2)):
        genesDB[tmpList[i*2][1:].rstrip()]=tmpList[i*2+1].rstrip();
    tmpList=[];
    return genesDB;

def curateDB(genesDB):
    curGenesDB={}
    scoreDB={}
    exclude=['rrlA', 'rrlB', 'rrlC', 'rrlD', 'rrlE', 'rrlF', 'rrlG', 'rrlH', 'rrfA', 'rrfB', 'rrfC', 'rrfD',
    'rrfE', 'rrfF', 'rrfG', 'rrfH', 'rrsA', 'rrsB', 'rrsC', 'rrsD', 'rrsE', 'rrsF', 'rrsG', 'rrsH',
    'rnpB','cysT','aspT','aspU','aspV','serV','serW','serX','serT','serU','glnV','glnX','glnU',
    'glnW','metV','metW','metY','metZ','metT','metU','asnT','asnU','asnV','asnW','proK','proM',
    'proL','lysQ','lysT','lysV','lysW','lysY','lysZ','thrT','thrV','thrW','thrU','pheU','pheV',
    'alaT','alaU','alaV','alaW','alaX','glyV','glyW','glyX','glyY','glyU','glyT','ileT','ileU',
    'ileV','leuP','leuQ','leuT','leuV','leuZ','leuX','leuU','leuW','hisR','argQ','argV','argY',
    'argZ','argX','argU','argW','trpT','valT','valU','valX','valY','valZ','valW','valV','gltT',
    'gltU','gltV','gltW','tyrT','tyrU','tyrV']
    for i in genesDB:
        if (not ((len(genesDB[i])<96) or (i in exclude) or len(genesDB[i])%3!=0)):
            curGenesDB[i]=genesDB[i];
    return curGenesDB;

def sortDB(genesDB):
    scoreDB={}
    sortGeneDB={}
    for i in genesDB:
        seqNew=str(genesDB[i])
        score=0
        codon=''
        for j in range(int(len(seqNew)/3)):
            codon=seqNew[3*j]+seqNew[3*j+1]+seqNew[3*j+2]
            score=score+codonDict[codon]
        score=score/(len(seqNew)/3)
        scoreDB[i]=score
    tmpSrt=sorted(scoreDB.items(), key = lambda kv:(kv[1], kv[0]));
    for i in tmpSrt:
        sortGeneDB[i[0]]=genesDB[i[0]]
    return sortGeneDB;

def Seq2Binary(seqIn):
    seqNew=seqIn;
    binSeq=[]
    a=[];
    for i in range(int(len(seqNew)/3)):
        codon=seqNew[3*i]+seqNew[3*i+1]+seqNew[3*i+2]
        a.append(codonDict[codon]);
        binSeq.append(a)
        a=[]
    return binSeq;

def Bin2Bytes(binSeqIn, maxGeneLen):
    a=[255]
    tmpBinSeq=binSeqIn
    for i in range(int(maxGeneLen)-int(len(binSeqIn))):
        tmpBinSeq.append(a)
    seqBytes=bytes();
    for i in tmpBinSeq:
        seqBytes=seqBytes+bytes(i);
    if len(seqBytes)==int(maxGeneLen):
        return seqBytes;
    else:
        print("error")
        exit(0)
        return seqBytes;

def dnaDCT(binSeqIn, useMed, dctType):
    binSeqList=[]
    for i in binSeqIn:
        for j in i:
            binSeqList.append(j)
    npSeqBin=np.array(binSeqList)
    npDCTseqBin=dct(npSeqBin, dctType)
    fingerprint=[]
    k=0
    for i in npDCTseqBin:
        if k<32:
            fingerprint.append(i)
        k=k+1
    medFP=statistics.median(fingerprint)
    dctHash=''
    for i in fingerprint:
        if useMed==0:
            if i > 0:
                dctHash=dctHash+'1'
            elif i<=0:
                dctHash=dctHash+'0'
        elif useMed==1:
            if i > medFP:
                dctHash=dctHash+'1'
            elif i<= medFP:
                dctHash=dctHash+'0'
    return dctHash;

def main():
    seqBytesM={}
    seqBinM={}
    seqByteFill={}
    sortSeqByteFill={}
    hashDB={}
    maxGeneLen=0
    fullPic=bytes();
    fName='/Users/trsheph/Documents/GitHub/SynBioStandardizer/SynColi3/MG1655_genes.txt'
    genesDB=sortDB(curateDB(loadGenomeData(fName)));
    k=1
    for i in genesDB:
        if len(genesDB[i])>maxGeneLen:
            maxGeneLen=len(genesDB[i]);
    maxGeneLen=int(maxGeneLen/3)
    for i in genesDB:
        seqBinM[i]=Seq2Binary(genesDB[i]);
        hashDB[i]=dnaDCT(seqBinM[i],1,1);
        fullPic=fullPic+Bin2Bytes(seqBinM[i], maxGeneLen);
        print(str(k)+'/'+str(len(genesDB)))
        k=k+1
    fOut=open('HashTable.txt', 'w')
    for i in hashDB:
        fOut.write(i+"\t"+hashDB[i]+"\n");
    fOut.close()
    seqImgC=Image.frombytes("L", (maxGeneLen,len(genesDB)), fullPic);
    seqImgC.save('SortByUsage.png')
    hammingWrite(hashDB)

if __name__ == "__main__":
    main();
