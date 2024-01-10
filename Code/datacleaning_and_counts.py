#数据清洗
##1 去除非atcg所在地方的三个字符
##2 去除完全包含的值
import pandas as pd
import re
def data_cleaning(loc):
    with open(loc, "r") as f:  # 读取文件
        lines = f.readlines()
    for n, line in enumerate(lines):  # 去除\n
        lines[n] = line.rstrip()
    all_seq = {}  # 储存序列
    i = 0  # 设置一个控制参数i,初始值为0
    for n, line in enumerate(lines):  # 遇到>时，i+1，i的值为1且这一行不包含>时，将这行写入line对应的值中
        if i == 1 and ">" not in line:
            all_seq[discription_line] = all_seq[discription_line] + lines[n]  # 写入，换行也继续写入
        if ">" in line:  # 当遇到序列描述行时，将该描述行对应的值设置为空字符串，并且将i+1，并且用discription_line储存描述行
            i += 1
            all_seq[line] = ""
            discription_line = line  # 储存>seq
        if i == 2:  # 当i等于2时，说明遇到了下一个序列描述行，因此将i重新设置为1
            i = 1  # 将i重新赋值为1
            continue  # 进入下一循环
    #all_seq储存了所有seq对应的密码子
    # print(len(all_seq))#原序列长度
    ##保证键不重复
    list_keys = list(all_seq.keys())  # 用列表储存所有的键
    all_keys_set = set(all_seq.keys())  # 将所有的键转换为集合
    # 检查集合的长度是否等于原始键的数量
    if len(all_keys_set) == len(list_keys):
        # print("所有键均不重复！")

        # 去除含有非ATCG的字符集合所在的密码子
        valid_chars = set("ATCG")  # 设置一个包含合法字符的集合
        remove_list = [] #储存需要移除的键
        for description, sequence in all_seq.items():
            invalid_chars = set(sequence) - valid_chars  # 计算非法字符
            if invalid_chars:
                # print(f"在{description}中发现非法字符: {', '.join(invalid_chars)}") #不合格的键值对描述
                remove_list.append(description)#储存不符合要求的键

        for key in remove_list:    #去除不符合要求的密码子对应的三个字符
            # print(all_seq[key])
            # print(len(all_seq[key]))
            value=""#value:用于记录新的密码子序列
            for i in range( 0, len(all_seq[key])-2 , 3):#每三个字符截取
                codon = all_seq[key][i:i+3]
                invalid_chars = set(codon) - valid_chars  # 计算非法字符
                if invalid_chars: #如果发现存在非法字符，则去除其所在的三个字符
                    # print(codon)
                    continue #不记录非法字符所在的三个密码子序列
                else:
                    value+=codon
            all_seq[key] = value#更新键值对
        ###all_seq:记录了所有去除非法字符的键值对

        removelist_genekey = []#用removelist_genekey储存需要删除的键
        genelist = [] #genelist储存从键中提取出的所有的gene
        repeatlist=[]
        ###去除gene相同时 完全包含密码子的情况

        ##先找出所有有重复gene的键值对
        for key in all_seq.keys():#遍历所有的键
            gene_match = re.search(r'gene:(\S+?)\s+gene_biotype:',key)#找出gene
            if gene_match:
                gene =gene_match.group(1)
                genelist.append(gene)
        # print(genelist)
        # list(all_seq.keys())列表记录所有的键
        # list(all_seq.values())列表记录所有的值

        ii=0#索引
        jj=0#标识，判断是否终止
        list0=[]
        for ii in range(len(genelist)):
            gene=genelist[ii]
            if genelist.count(gene) !=1:#如果多次出现该基因
                jj = 0
                repeatlist.append(ii)
            else:
                jj = 1
            if jj == 1 and len(repeatlist)>0:
                list0.append(repeatlist)
                repeatlist=[]
            ii+=1

        # list0储存了所有重复的gene序列所在的序号

        for i in list0:
            list_compare = []#list_compare 储存所有gene相同时的序列
            for j in i:
                list_compare.append(list(all_seq.values())[j])
            unique_items=[]
            for item in list_compare:
                is_unique = all(item not in other_item for other_item in list_compare if other_item != item)
                if is_unique:
                    unique_items.append(item)
            result = [item for item in list_compare if item not in unique_items]
            ##找出需要删除的key
            for k, v in all_seq.items():
                if v in result:
                    removelist_genekey.append(k)
        #删除重复的key
        for key in removelist_genekey:
            if key in all_seq:
                del all_seq[key]
        # print(len(all_seq))

    return all_seq


def calculate_gene_counts(gene_dict, codon_table_path, output_result_path):
    # gene_dic 是数据清洗后得到的结果，
    # codon_table_path 是密码子氨基酸对应关系表格，放在相对路径下，
    # output_result_path是输出结果保存的路径。可直接输入文件名保存在相对路径下。

    # 读取密码子和蛋白质对应关系表
    codon_table_df = pd.read_csv(codon_table_path)

    # 初始化最终的df,增加了两列codon_counts和protein_counts，并赋值为0
    result_df = codon_table_df.copy()
    result_df['codon_counts'] = 0
    result_df['protein_counts'] = 0
    # print(result_df)

    # 初始化总计数的字典
    total_codon_counts = {}
    total_protein_counts = {}

    # 遍历基因字典
    for gene, dna_sequence in gene_dict.items():
        protein_sequence = ''
        codon_counts = {}
        protein_counts = {}

        # 遍历DNA序列并翻译
        for i in range(0, len(dna_sequence) - len(dna_sequence) % 3, 3):
            codon = dna_sequence[i:i + 3]
            #  在codon_table_df中查找密码子对应的蛋白质。
            protein = codon_table_df.loc[codon_table_df['codon'] == codon, 'protein1'].values[0]
            # print(codon_table_df.loc[codon_table_df['codon'] == codon, 'protein1'])  # 此句代码返回的是对应密码子所在的行和氨基酸。
            protein_sequence += protein

            # 更新密码子计数
            codon_counts[codon] = codon_counts.get(codon, 0) + 1  # 如果字典没有codon，则初始化为0，再+1.如果键存在就得到相应的值，再+1

            # 更新氨基酸计数
            protein_counts[protein] = protein_counts.get(protein, 0) + 1
        print(total_codon_counts)
        print(total_protein_counts)
        # 将当前基因的计数累加到总计数中，将一条DNA的结果加到之前的结果上。
        for codon, count in codon_counts.items():
            total_codon_counts[codon] = total_codon_counts.get(codon, 0) + count

        for protein, count in protein_counts.items():
            total_protein_counts[protein] = total_protein_counts.get(protein, 0) + count

    # 更新总计数到result_df中，根据DataFrame中的‘codon’列和‘protein1’列将字典和DataFrame匹配起来。并把NA填成0
    result_df['codon_counts'] = result_df['codon'].map(total_codon_counts).fillna(0)
    result_df['protein_counts'] = result_df['protein1'].map(total_protein_counts).fillna(0)

    # 保存结果
    result_df.to_csv(output_result_path, index=False)

    return result_df

loc_Arabidopsis = "../Data/Fasta/Arabidopsis_thaliana.TAIR10.cds.all.fa"
loc_C_elegans = "../Data/Fasta/Caenorhabditis_elegans.WBcel235.cds.all.fa"
loc_Drosophila = "../Data/Fasta/Drosophila_melanogaster.BDGP6.46.cds.all.fa"
loc_human = "../Data/Fasta/Homo_sapiens.GRCh38.cds.all.fa"
loc_mouse = "../Data/Fasta/Mus_musculus.GRCm39.cds.all.fa"
loc_yeast = "../Data/Fasta/Saccharomyces_cerevisiae.R64-1-1.cds.all.fa"


Arabidopsis = data_cleaning(loc_Arabidopsis) #拟南芥字典
C_elegants = data_cleaning(loc_C_elegans) #线虫字典
Drosophila = data_cleaning(loc_Drosophila) #果蝇字典
Human = data_cleaning(loc_human) #人字典
Mouse = data_cleaning(loc_mouse) #小鼠字典
Yeast = data_cleaning(loc_yeast) #酵母字典

# print(Arabidopsis)

codon_table_path = 'codon_table.csv'
Human_result_df = calculate_gene_counts(Human, codon_table_path, '../Result/Counts/Human_counts.csv')
Arabidopsis_result_df = calculate_gene_counts(Human, codon_table_path, '../Result/Counts/Arabidopsis_counts.csv')
C_elegants_result_df = calculate_gene_counts(Human, codon_table_path, '../Result/Counts/C_elegants_counts.csv')
Drosophila_result_df = calculate_gene_counts(Human, codon_table_path, '../Result/Counts/Drosophila_counts.csv')
Mouse_result_df = calculate_gene_counts(Human, codon_table_path, '../Result/Counts/Mouse_counts.csv')
Yeast_result_df = calculate_gene_counts(Human, codon_table_path, '../Result/Counts/Yeast_counts.csv')
