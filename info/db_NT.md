## NCBI Nucleotide database (nt) wihtout viruses
The following are the steps to construct non-viral nt database.

**Steps bellow assume that you have set `VirID_DB_PATH` path**

Tools Requirements: [`BLAST+`](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) and [`seqkit`](https://github.com/shenwei356/seqkit/releases)

---


## **Method 1**: Download and configure from scratch an nt database that does not contain virus sequences.


### Step 1.1 Download the pre-indexed NCBI BLAST non-viral nt database

```shell
cd $VirID_DB_PATH

mkdir $VirID_DB_PATH/NT

# Eukaryota nt
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_euk-nucl-metadata.json

# nt_euk-nucl.filesList: URL for downloading eukaryotic data
grep "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt_euk" nt_euk-nucl-metadata.json |sed 's/[[:space:]]*["|,]//g' > nt_euk-nucl.filesList

# downloading, use `aspera` tool rather than `wget` might be faster.
cat nt_euk-nucl.filesList |
while read line;
do
wget -c "$line";
wget -c "$line".md5;
done
```


```shell
# Prokaryota (bacteria and archaea) nt
wget https://ftp.ncbi.nlm.nih.gov/blast/db/nt_prok-nucl-metadata.json

# nt_prok-nucl.filesList: URL for downloading Prokaryota data
grep "ftp://ftp.ncbi.nlm.nih.gov/blast/db/nt_prok" nt_prok-nucl-metadata.json |sed 's/[[:space:]]*["|,]//g' > nt_prok-nucl.filesList

# Use `aspera` tool rather than `wget` might be faster.
cat nt_prok-nucl.filesList |
while read line;
do
wget -c "$line";
wget -c "$line".md5;
done
```


`aspera` example :
```shell
for i in {00..108};
do 
ascp -QT -i ~/.aspera/connect/etc/asperaweb_id_dsa.openssh -l 100M -k 1 -T anonftp@ftp.ncbi.nlm.nih.gov:/blast/db/nt_euk."$i".tar.gz ./
done
```
---
### Step 1.2 Check integrity


```shell
md5sum -c *.md5
```
If the all files are ok, then next step, or you need to download again.


---
### Step 1.3 Decompress


```shell
cat nt_prok-nucl.filesList nt_euk-nucl.filesList |
while read line;
do
tar -zxvf "$line";
done

```
### Step 1.4 Aggregate existing nt_prok & nt_euk database


```shell

blastdb_aliastool -dblist "nt_prok nt_euk" -dbtype nucl -out $VirID_DB_PATH/NT/nt -title "NT database with only Eukaryota & Prokaryota (bacteria and archaea)" 

```

Method 1 finished here.


---

## **Method 2**: The complete fasta sequence of the nt database already exists.

***If you haven't downloaded the fasta sequence file of the nt database yet, we strongly recommend using method 1.***


Suppose the FASTA file name of your nt database is: `nt.fasta`, and it has been moved to `$VirID_DB_PATH/NT/nt.fasta`.


```shell
cd $VirID_DB_PATH/NT
```


---
### Step 2.1. Prepare the viral sequences

#### Setp 2.1.1 Prepare the viral sequences accession number listed in **NCBI Virus**.
- We have pre-collected [all NCBI virus nucleic acid sequence accession numbers](./NCBI_Virus_AccessionList_20240617.tar.gz) available in the NCBI Virus website.(as of 17 June 2024) (Another virus accession list download [link](https://drive.google.com/file/d/1Rx_CGeg0yPOAfMJwVCxMVFFNgBPzOo3r/view?usp=drive_link))


```shell
# Download manually from the link above, and move it into $VirID_DB_PATH/NT

tar -jxvf NCBI_Virus_AccessionList_20240617.tar.gz

cat NCBI_Virus_*_AccessionList_20240617.tsv > NCBI_Virus_All_Nucleotides_Accession_List_20240617.tsv
```


#### Setp 2.1.2 Prepare the viral sequences accession number listed in **NCBI nt virus database**.
- We also pre-collected [nt virus accession number list](https://drive.google.com/file/d/1X2R4a0pDNIOiNSRwPDa6x_Wy95eYHo4E/view?usp=sharing).(as of 11 July 2024) 

```shell
# Download manually from the link above, and move it into $VirID_DB_PATH/NT

tar -jxvf nt_viruses_acc.tsv.tar.gz
```

#### Setp 2.1.3 Prepare the artificial and other sequences accession number list.
- Download the pre-collected list [here](https://drive.google.com/file/d/1s9T-vIn3zQLndEah20WnDJysRz6Mtn3E/view?usp=sharing).(as of 11 July 2024) and move it into `$VirID_DB_PATH/NT`


Then:


```shell
cat NCBI_Virus_All_Nucleotides_Accession_List_20240617.tsv nt_viruses_acc.tsv nt_others_acc.tsv > VirID_nt_exclude_acc.tsv
```


---


### Step 2.2 Remove the viral sequences

```shell

seqkit grep -f VirID_nt_exclude_acc.tsv -v nt.fasta -o nt_WO_Virus.fasta

```

---

### Step 2.3 Build the index of nt_WO_Virus.fasta for `blastn`.

```shell

makeblastdb -in nt_WO_Virus.fasta -dbtype nucl  -hash_index -out $VirID_DB_PATH/NT/nt -parse_seqids

```

Method 2 finished here.

