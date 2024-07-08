## NCBI Non-Redundant Protein Database (NR)

We use diamond to perform protein sequence alignment.

There are two methods to construct the necessary NR database for use with diamond

Here we provide two ways to download the NR database, but there are also other ways to download it.

Files are large and can take a lot of time depending on the network.

**Note:The final index path for both methods needs to be: `"$VirID_DB_PATH"/NR/nr`**

### 1. From the sequences

  - **1.1 Download the sequence file.**
    ```shell
      wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz

      wget -c https://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz.md5
    ```

  - **1.2 Check for the file integrity.**
    ```shell
      md5sum -c nr.gz.md5
    ```
      If the result is ok, then next step

  - **1.3 Using ***[diamond](https://github.com/bbuchfink/diamond)*** to build the index.**
    ```shell
      diamond makedb --in nr.gz -d "$VirID_DB_PATH"/NR/nr
    ```

### 2. From the pre-builded BLAST index

  - **2.1 Download the index file.**

     These indexes are constructed from multiple files, currently numbered from 00 to 97, and possibly more in the future.(from `nr.00.tar.gz` to `nr.97.tar.gz`)

      ```shell
      cd VirID_DB_PATH/NR

      for i in {00..97}; 
        do
          base_url="https://ftp.ncbi.nlm.nih.gov/blast/db/nr.$i"
          tar_url="$base_url.tar.gz"
          md5_url="$tar_url.md5"

          wget $tar_url
          wget $md5_url
        done

      ```


  - **2.2 Check for the file integrity.**
    ```shell
    for i in {00..97}; 
      do
        md5_path="nr.$i.tar.gz.md5"
        md5sum -c "$md5_path"
      done
    ```
    If the results is ok, then next step

  - **2.3 Unzip the files.**
    ```shell
    for i in {00..97}; 
      do
        index_path="nr.$i.tar.gz"
        tar -zxvf "$index_path" -C "$VirID_DB_PATH"/NR
      done
    ```

  - **2.4 Rebuild the indexes.**
    
    **Note:**
  To use the *prepdb* command, you've been advised to install [diamond]((https://github.com/bbuchfink/diamond)) from its GitHub repository instead of using the Conda installation.
      
      ```shell
      diamond prepdb -d VirID_DB_PATH/NR/nr
      ```


