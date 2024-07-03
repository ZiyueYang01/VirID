import dendropy
import os
import logging
import pandas as pd
from Bio import SeqIO
from os import system
from VirID.external.blast import Blastp
from VirID.config.config import ICTV_TAXON_CSV

class Tree_taxon(object):
    
    def __init__(self,input_tree,ref_tsv,clustr_name,item,msa,out_dir):
        self.input_tree = input_tree
        self.ref_tsv = ref_tsv
        self.clustr_name = clustr_name
        self.item = item
        self.ref_id = list(pd.read_csv(self.ref_tsv)['Cluster_sseqid'])
        self.csvPD = pd.read_csv(self.ref_tsv,sep=',',header=0)
        self.msa = msa
        self.out_dir = out_dir
        self.logger = logging.getLogger('timestamp')



#这里我想构建一个字典，保存一下，不能每次都进行blastp吧，太浪费时间了
    def find_ref_class(self,node_ID):
        self.csvPD = pd.read_csv(self.ref_tsv,sep=',',header=0)
        for i in range(len(self.csvPD)):
            if str(self.csvPD['Cluster_sseqid'][i]) == node_ID.strip("'"):
                return(self.csvPD[self.item][i])
        return None


    def read_tree(self):
        tree = dendropy.Tree.get_from_path(self.input_tree,
                                        schema='newick',
                                        rooting='force-unrooted',
                                        preserve_underscores=True)
        tree.reroot_at_midpoint()
        tree.ladderize(ascending=False)
        return tree

    def cal_aai(self,node,sister):
        supergroup = self.out_dir.split("/")[-1]
        out_fas = os.path.join(self.out_dir,supergroup)+"_all.fasta"

        seq_trim = SeqIO.to_dict(SeqIO.parse(out_fas, "fasta"))
        query_str = str(seq_trim[str(node.taxon).strip("'")].seq)
        ref_str = str(seq_trim[str(sister.taxon).strip("'")].seq)

        query_str_file = f'{self.out_dir}/blastp.fas'
        ref_str_file = f'{self.out_dir}/blastp_ref.fas'
        out_tsv = f'{self.out_dir}/blastp.tsv'

        with open(query_str_file, 'w') as fw:
            fw.write(f">query"+ '\n')
            fw.write(str(query_str) + '\n')
        with open(ref_str_file, 'w') as fw:
            fw.write(f">ref"+ '\n')
            fw.write(str(ref_str) + '\n')

        db_path = f'{self.out_dir}/blastp'
        system(f'makeblastdb \
                -dbtype prot \
                -in {ref_str_file} \
                -input_type fasta \
                -parse_seqids \
                -out {db_path} ')
        
        aa_blastp = Blastp().run(query_str_file,"blastp",out_tsv,db_path)

        size = os.path.getsize(aa_blastp)
        if size <1 :
            return 0

        aa_table = pd.read_csv(aa_blastp,header=None,encoding = "utf-8",sep = '\t', 
                                names =['qaccver','saccver','pident','length',"mismatch","gapopen",
                                        'qstart','qend','sstart','send','evalue','bitscore','qlen','slen'])
        
        return aa_table['pident'][0]

    def kids_has_ref(self,node):
        if str(node.taxon).strip("'") in self.ref_id:
            return True
        res = []
        for child in node.child_node_iter(): 
            sister_node = child.sister_nodes()[0]
            if (str(child.taxon).strip("'") in self.ref_id) and ((str(sister_node.taxon).strip("'") in self.ref_id) is False and sister_node.is_leaf()):
                child.label = self.find_ref_class(str(child.taxon))

                if self.item != 'Human':
                    sister_node.label = child.label
                else:
                    identity = self.cal_aai(child,sister_node)
                    if int(identity) < 80:
                        sister_node.label = self.clustr_name

                    else:
                        sister_node.label = child.label
                node.label = sister_node.label
                self.par_node_class(node)
                return True
            if child.is_internal():
                res.append(self.kids_has_ref(child))
            if str(child.taxon).strip("'") in self.ref_id:
                res.append(True)
                return True
            else:
                res.append(False)
        res = list(set(res))
        if(True in res):
            return True
        else:
            return False


    def sister_node(self,node):
        sister = node.sister_nodes()[0]
        if str(node.label) != 'None':
            return node
        elif (str(sister.taxon).strip("'")  in self.ref_id):
            return sister
        elif sister.is_internal(): 
            if self.kids_has_ref(sister) is False:
                sister = self.sister_node(sister.parent_node)
            else:
                return sister
        else:
            sister = self.sister_node(sister.parent_node)
        return sister


    def long_str(self,li): 
        res=''    
        num_rows = len(li)
        for i in li:
            if len(i)>11:
                del i[11:]
        num_cols = min(len(row) for row in li)
        result = []
        for col in range(num_cols):
            common_element = li[0][col]
            for row in range(1, num_rows):
                if li[row][col] != common_element:
                    return "|".join(result)
            result.append(common_element)
        if len(result) < 11:
            for i in range(11-len(result)):
                result.append("")
        res =  "|".join(result)
        return res


    def decorate_inter_node(self,node):
        if str(node.taxon).strip("'") in self.ref_id:
            node.label = self.find_ref_class(str(node.taxon))
            return node.label
        if node.is_leaf() or str(node.label) != 'None':
            return node.label
        if str(node.label) ==  'None':
            name_set = []
            for child in node.child_node_iter():
                if child.is_internal():
                    child.label = self.decorate_inter_node(child)
                if (str(child.taxon).strip("'") in self.ref_id) is False and child.is_leaf():
                    child_sister = child.sister_nodes()[0]
                    new_node = self.sister_node(child_sister)
                    child.label = self.decorate_inter_node(new_node)
                if self.item == 'Taxon_set':
                    name_set.append(str(child.label).split("|"))
                else:
                    name_set.append(str(child.label))
            if self.item == 'Taxon_set':
                max_str = self.long_str(name_set)   
                return max_str
            else:
                name_set = list(set(name_set))
                if(len(name_set) == 1):
                    return name_set[0]
                else:
                    return self.clustr_name
        return node.label


    def par_node_class(self,node):
        sister = node.sister_nodes()[0]
        if str(node.label) != 'None' and str(sister.label) != 'None':
            par = node.parent_node
            if node.label == sister.label:
                par.label = node.label
            else:
                if self.item == 'Taxon_set':
                    par.label = self.long_str([str(node.label).split("|"),str(sister.label).split("|")])
                else:
                    node_label = str(node.label).split("_")[0] if ('_' in str(node.label)) else str(node.label)
                    sister_label = str(sister.label).split("_")[0] if ('_' in str(sister.label)) else str(sister.label)
                    par.label = node_label  if (node_label == sister_label) else self.clustr_name


    def decorate(self,tree):
        for node in tree.internal_nodes():
            node.label = None

        for node in tree.leaf_node_iter():
            if str(node.label) ==  'None':
                sister = node.sister_nodes()[0]
                if str(node.taxon).strip("'") in self.ref_id: 
                    taxon = self.find_ref_class(str(node.taxon))
                    if taxon is not None:
                        node.label = taxon
                    else:
                        node.label = self.clustr_name
                    if sister.is_internal(): 
                        sister.label = self.decorate_inter_node(sister)
                        self.par_node_class(node)
                    elif str(sister.taxon).strip("'") in self.ref_id:
                        sister.label = self.find_ref_class(str(sister.taxon))     
                        self.par_node_class(node)
                elif str(sister.taxon).strip("'") in self.ref_id:
                    taxon = self.find_ref_class(str(sister.taxon))
                    if taxon is not None:
                        sister.label = taxon
                    else:
                        sister.label = self.clustr_name
                    if self.item != 'Human':
                        node.label = sister.label
                    else:
                        identity = self.cal_aai(node,sister)
                        if int(identity) < 80:
                            node.label = self.clustr_name
                        else:
                            node.label = sister.label
                    self.par_node_class(node)
                elif str(node.taxon).strip("'") not in self.ref_id and sister.is_internal():
                    if self.kids_has_ref(sister):
                        new_node = self.sister_node(node)
                        node.label = self.decorate_inter_node(new_node)
                        self.par_node_class(node)

        for node in tree.leaf_node_iter():
            if str(node.label).strip("'") ==  'None':
                new_node = self.sister_node(node)
                node.label = self.decorate_inter_node(new_node)
                if str(node.label) ==  'None':
                    node.label = self.clustr_name
                self.par_node_class(node)



    def classify(self,tree,query_list):
        query_dic = {}
        taxon_dic = {}
        all_group_csv = pd.read_csv(ICTV_TAXON_CSV,sep=',',header=0)
        all_group_pd = pd.DataFrame(all_group_csv)
        ictv = ['Realm','Kingdom','Phylum','Subphylum','Class','Order','Suborder','Family','Subfamily','Genus','Species']
        for node in query_list:
            node_taxon = str(node.taxon).strip("'")
            query_node = tree.find_node_with_taxon_label(node_taxon)
            query_dic[node_taxon] = query_node.label

            sister = query_node.sister_nodes()[0]
            if self.item == 'Taxon_set':
                if str(sister.taxon).strip("'") in self.ref_id: 
                    itcv_taxon = []
                    query_dic[node_taxon] = str(self.csvPD[self.csvPD['Cluster_sseqid']==str(sister.taxon).strip("'")]["Cluster_taxon"].values).strip("[").strip("]").strip("'")
                    for a in ictv:
                        itcv_taxon.append(str(self.csvPD[self.csvPD['Cluster_sseqid']==str(sister.taxon).strip("'")][a].values).strip("[").strip("]").strip("'"))
                    itcv_taxon.append(str(sister.taxon).strip("'"))
                    taxon_dic[node_taxon] = itcv_taxon
                else:
                    if str(node.label) ==  '':
                        itcv_taxon = []
                        for a in ictv:
                            itcv_taxon.append(str(all_group_pd[all_group_pd['TAXON']==self.clustr_name][a].values).strip("[").strip("]").strip("'"))
                        itcv_taxon.append('')
                        taxon_dic[node_taxon] = itcv_taxon
                    else:
                        itcv_taxon = str(query_node.label).split("|")

                        closest_string = ""
                        all_group_pd_part = all_group_pd[all_group_pd['RdRP_super_group']== self.clustr_name]
                        all_group_pd_group = all_group_pd_part['Taxon_set'].values
                        for ref_ID in all_group_pd_group:
                            res_list = self.long_str([str(ref_ID).split("|"),str(query_node.label).split("|")]).split("|")
                            new_list = [x for x in res_list if x != "|"]
                            closest_string = "|".join(new_list) if len(new_list) > len(closest_string.split("|"))  else closest_string

                        if len(str(closest_string).split("|")) < 11:
                            for i in range(11-len(str(closest_string).split("|"))):
                                closest_string += "|"

                        query_dic[node_taxon] = str(all_group_pd_part[all_group_pd_part['Taxon_set']==str(closest_string)]['TAXON'].values).strip("[").strip("]").strip("'")
                        for i in range(12-len(itcv_taxon)):
                            itcv_taxon.append("")
                        taxon_dic[node_taxon] = itcv_taxon
        return query_dic,taxon_dic


    def run(self):
        query_list = []
        tree = self.read_tree()

        ref_ID = self.ref_tsv
        ref_id = list(pd.read_csv(ref_ID)['Cluster_sseqid'])

        for node in tree.leaf_node_iter():
            if (str(node.taxon).strip("'") in ref_id) is False:
                query_list.append(node)

        self.decorate(tree)
        out_label_tree = os.path.join(self.out_dir,self.clustr_name)+"_query_label_file.txt"
        
        if(os.path.isfile(out_label_tree)):
            os.remove(out_label_tree)

        tree.write(path=out_label_tree,schema="nexus")

        query_dic,taxon_dic = self.classify(tree,query_list)
        query_dic_csv = pd.DataFrame(list(query_dic.items()))

        if self.item == 'Taxon_set':
            taxon_dic_csv = pd.DataFrame.from_dict(taxon_dic,orient = 'index')
            taxon_dic_csv = taxon_dic_csv.reset_index(drop=False)
            taxon_dic_csv.columns = ['qseqid','Realm','Kingdom','Phylum','Subphylum','Class','Order','Suborder','Family','Subfamily','Genus','Species','RdRP_Species']
            return query_dic_csv,taxon_dic_csv
        else:
            return query_dic_csv
