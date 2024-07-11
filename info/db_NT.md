## NCBI Nucleotide database (nt) wihtout viruses
The following are the steps to exclude virus sequences from the nt database from the sequence file in the `fasta` format 

**The following steps usually assumes that all files are in the same directory. When referring to actual use, please pay attention to the file path**


### Step 1. Download the NT database in FASTA format
[NT database in FASTA format](https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz)  over 378 GB
  - **1.1 Download the sequence file.**
    ```shell
      wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz

      wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nt.gz.md5
    ```
    Files are large and can take a lot of time depending on the network.
    
    **Other download methods may be faster,`aspera` et al.**

  - **1.2 Check for the file integrity.**
    ```shell
      md5sum -c nt.gz.md5
    ```
      If the result is ok, then next step
    
  - **1.3 Decompress.**
    ```shell
      gunzip nt.gz
    ```
    After decompression, you will get `nt.fasta`
    

### Step 2. Remove the viral sequences listed in [NCBI Virus](https://www.ncbi.nlm.nih.gov/labs/virus/vssi/#/).
We have pre-collected all virus nucleic acid [sequence IDs](./NCBI_Virus_All_Nucleotides_Accession_List_20240617.tar.gz) (as of 17 June 2024) available in the NCBI Virus website.

  - **2.1 Download and decompress the virus nucleic acid sequence IDs list.**
    ```shell
      wget -c https://github.com/ZiyueYang01/VirID/blob/main/info/NCBI_Virus_All_Nucleotides_Accession_List_20240617.tar.gz
    
      tar -zxvf NCBI_Virus_All_Nucleotides_Accession_List_20240617.tar.gz
    
      cat NCBI_Virus_*_AccessionList_20240617.tsv > NCBI_Virus_All_Nucleotides_Accession_List_20240617.tsv
    ```
  - **2.2 Download and decompress the virus nucleic acid sequence IDs list.**
    ```shell
      seqkit grep -f NCBI_Virus_All_Nucleotides_Accession_List_20240617.tsv -v nt.fasta -o nt_noViruses_0.fas
    ```


### Step 3.Remove the viral sequences based on their taxonomic IDs (taxids)


  - **3.1 Download these files.**

    - [nucl_gb.accession2taxid](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz)  2.3 GB

    - [taxdump](https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz) 62.7 MB
   
    Example:
    ```shell
      wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/accession2taxid/nucl_gb.accession2taxid.gz
    
      wget -c https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz

    ```
    Files are large and can take a lot of time depending on the network.
    
    **Other download methods may be faster,`aspera` et al.**


#### Assuming you have downloaded them all to your Linux directory now: `/this/path/`

```shell

cd /this/path/

gunzip nucl_gb.accession2taxid.gz
#nucl_gb.accession2taxid

tar -zxvf taxdump.tar.gz
#a lot of files

seqkit seq -n -i NT_noViruses_0.fas > NT_noViruses_0.id

sed -i '/^$/d' NT_noViruses_0.id

grep -F -f NT_noViruses_0.id nucl_gb.accession2taxid > NT_noViruses_0.id2taxid

sed -i '/^$/d' NT_noViruses_0.id2taxid


cut -f3 NT_noViruses_0.id2taxid \
| sort -u \
| taxonkit --data-dir /this/path/ lineage \
| awk '$2>0' \
| taxonkit --data-dir /this/path/ reformat -f "{k}\t{kID}\t{p}\t{pID}\t{c}\t{cID}\t{o}\t{oID}\t{f}\t{fID}\t{g}\t{gID}\t{s}\t{sID}" -F \
| cut -f1,2,4,6,8,10,12,14 \
> NT_noViruses_0.tax

awk -F'\t' '$2 != 10239' NT_noViruses_0.tax > NT_noViruses_1.tax

```
