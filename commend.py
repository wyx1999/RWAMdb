import json
import os
import django
from django.db.models import Count, Q
import pandas as pd
from django.db.models import F, Func, Value
from django.db.models.functions import Upper, Lower
from django.db import transaction
# 设置 Django 的环境变量
os.environ.setdefault("DJANGO_SETTINGS_MODULE", "lab.settings")
django.setup()
import math
from django.db.models.functions import Concat
# 现在可以使用 Django 的 ORM 和其他功能
from db.models import MAG, Gene ,Site ,Sample, Pathway

def creat_database():
    """
    Site.objects.all().delete()
    Sample.objects.all().delete()

    sample_df = pd.read_excel('db/static/db/sample.xlsx', header=0)
    sample_df = sample_df.values.tolist()
    site_df = pd.read_excel('db/static/db/site.xlsx', header=0)
    site_df = site_df.values.tolist()


    for i in site_df:
        Site.objects.create(name=i[0], longitude=i[1], latitude=i[2], altitude=i[3])


    for i in sample_df:
        Sample.objects.create(name=i[0], site=Site.objects.get(name=i[1]), water=i[2], TN=i[3], TP=i[4], TC=i[5], TOC=i[6], TS=i[7],MC=i[8], NO3=i[9], NH4=i[10], pH=i[11], ome= True if i[12]==1 else False)
    Bin.objects.all().delete()
    count = 0
    for i in bin_df:
        print(count)
        count += 1

        if str(i[5]) == 'nan':
            continue
        Bin.objects.create(name=i[4], sample=Sample.objects.get(name=i[1] + str(i[2])),
                           Kingdom=i[6].split('__')[1], Phylum=i[7].split('__')[1], Class=i[8].split('__')[1],
                           Order=i[9].split('__')[1], Family=i[10].split('__')[1], Genus=i[11].split('__')[1],
                           Species=i[12].split('__')[1],
                           relative_abundance=i[13], length=i[14], RPKM=i[15], Completeness=i[16], Contamination=i[17],
                           GC=i[18], N50=i[19])

    bin_df = pd.read_excel('db/static/db/bin.xlsx', header=0)
    bin_df = bin_df.values.tolist()

    site_df = pd.read_excel('db/static/db/site.xlsx', header=0)
    site_df = site_df.values.tolist()

    for i in site_df:
        Site.objects.create(name=i[0], longitude=i[1], latitude=i[2], altitude=i[3])

    t = set()
    for i in bin_df:
        t.add(i[1])
    print(len(t))



    bins = Bin.objects.all()
    data = []
    for bin in bins:
        a = [bin.name, bin.Taxonomy]
        data.append(a)
    df = pd.DataFrame(data, columns=['bin', 'taxonomy'])

    df.to_excel('db/static/db/taxonomy.xlsx', index=False)
    """
    # 打印当前工作目录（即程序运行时的目录）
    print("当前工作目录:", os.getcwd())


    # 修剪 Gene 表中 name 字段两端的空格
def batch_trim_gene_names(batch_size=10000):
    from django.db import connection
    total = Gene.objects.count()
    print(f"总记录数: {total}")

    for offset in range(0, total, batch_size):
        with connection.cursor() as cursor:
            cursor.execute(
                f"""
                UPDATE db_gene
                SET name = TRIM(name)
                WHERE rowid IN (
                    SELECT rowid FROM db_gene LIMIT {batch_size} OFFSET {offset}
                );
                """
            )
        print(f"已更新: {min(offset + batch_size, total)} / {total}")




def add_genes():
    with transaction.atomic():
        Gene.objects.all().delete()
    bins = MAG.objects.all()
    BATCH_SIZE = 10000
    for i, bin in enumerate(bins):
        count = 0
        print(bin.old)
        gene_data = []
        with open(f'F:/database file/annotation/{bin.old}/{bin.old}.emapper.genepred.fasta', 'r') as f:
            name = None
            for line in f:
                if line.startswith('>'):

                    # 如果之前的基因存在，插入 gene_data
                    if name:
                        count = count + 1
                        gene_data.append(
                            Gene(
                                MAG=bin,
                                id = bin.id+'_G'+str(count).zfill(5),
                                name=name,
                                length=sequence_length
                            )
                        )

                    name = line.split('>')[-1].strip()
                    sequence_length = 0  # 重置序列长度
                else:
                    # 只计算这一行的序列长度
                    sequence_length = len(line.strip())

            # 最后一行数据插入
            if name and sequence_length >= 0:
                count = count + 1
                gene_data.append(
                    Gene(
                        MAG=bin,
                        id=bin.id + '_G' + str(count).zfill(5),
                        name=name,
                        length=sequence_length
                    )
                )

        if gene_data:
            print(f"Inserting {len(gene_data)} genes for bin {bin.id}")
            with transaction.atomic():
                for i in range(0, len(gene_data), BATCH_SIZE):
                    Gene.objects.bulk_create(gene_data[i:i + BATCH_SIZE])


def add_pathway():
    c = 0
    files = ['M00165','M00173','M00374','M00375','M00376','M00377']
    dic = {'M00165': 'Reductive pentose phosphate cycle (Calvin cycle)',   # 19 11
           'M00173': 'Reductive citrate cycle (Arnon-Buchanan cycle)',
           'M00374': 'Dicarboxylate-hydroxybutyrate cycle',
           'M00375': 'Hydroxypropionate-hydroxybutylate cycle',
           'M00376': '3-Hydroxypropionate bi-cycle',
           'M00377': 'Reductive acetyl-CoA pathway (Wood-Ljungdahl pathway)'}
    gene_all = {'M00165': 19,  # 19 11
           'M00173': 40,
           'M00374': 24,
           'M00375': 19,
           'M00376': 21,
           'M00377': 7}
    step_all = {'M00165': 11,  # 19 11
                'M00173': 10,
                'M00374': 13,
                'M00375': 14,
                'M00376': 13,
                'M00377': 16}
    for file in files:
        df = pd.read_excel(f'db/static/db/pathway/{file}.xlsx')

        # 遍历每一行，将数据存入数据库
        for index, row in df.iterrows():
            # 获取每一行的数据

            bin_obj = MAG.objects.get(old=row['bin'])


            genes_number = row['genes_number']
            steps_number = row['steps_number']
            genes_percentage = row['genes_percentage']
            steps_percentage = row['steps_percentage']
            genes = row['genes']


            pathway = Pathway(
                bin=bin_obj,
                name=dic[file],
                module_ID=file,
                genes_number=genes_number,
                steps_number=steps_number,
                genes_percentage=genes_percentage,
                steps_percentage=steps_percentage,
                total_genes = gene_all[file],
                total_steps = step_all[file],
                genes=genes
            )
            pathway.save()
            c = c + 1
            print(c)

def patyway_update():
    Pathway.objects.filter(module_ID="M00165").update(
        total_genes=19,
        total_steps=11
    )
    Pathway.objects.filter(module_ID="M00173").update(
        total_genes=40,
        total_steps=10
    )
    Pathway.objects.filter(module_ID="M00374").update(
        total_genes=24,
        total_steps=13
    )
    Pathway.objects.filter(module_ID="M00375").update(
        total_genes=19,
        total_steps=14
    )
    Pathway.objects.filter(module_ID="M00376").update(
        total_genes=21,
        total_steps=13
    )

    Pathway.objects.filter(module_ID="M00377").update(
        total_genes=7,
        total_steps=16
    )
def add_site():
    df = pd.read_excel(f'db/static/db/site.xlsx')
    print(df)
    for _, row in df.iterrows():
        site_data = {
            "id": row["id"],
            "name": row["name"],
            "longitude": row.get("longitude", None),
            "latitude": row.get("latitude", None),
            "altitude": row.get("altitude", None),
            "temperature": row.get("temperature", None),
            "precipitation": row.get("precipitation", None),
            "temperature_2022": row.get("temperature_2022", None),
            "precipitation_2022": row.get("precipitation_2022", None),
            "type": row.get("type", None),
            "location": row["location"],
            "introduction": row.get("introduction", None),
            "chinese": row.get("chinese", None),
            "old": row.get("old", None),
        }

        # 使用 get_or_create 避免重复创建
        Site.objects.update_or_create(
            id=site_data["id"], defaults=site_data
        )

def add_sample():
    df = pd.read_excel(f'db/static/db/sample.xlsx')
    print(df)
    for index, row in df.iterrows():
        try:
            # 获取关联的 Site 对象
            site = Site.objects.get(id=row['site'])  # `site` 是 Excel 中的列名，确保其与数据库一致

            # 创建 Sample 对象
            sample = Sample(
                id=row['id'],
                site=site,
                TC=row.get('TC'),
                TS=row.get('TS'),
                TN=row.get('TN'),
                TP=row.get('TP'),
                SOC=row.get('SOC'),
                C_N=row.get('C_N'),
                IC=row.get('IC'),
                MBC=row.get('MBC'),
                pH=row.get('pH'),
                water=row.get('water'),
                NO3=row.get('NO3'),
                NH4=row.get('NH4'),
                date=row.get('date'),
                person=row.get('person'),
                old=row.get('old'),
            )
            sample.save()
            print(f"样本 {row['id']} 已成功保存！")
        except Site.DoesNotExist:
            print(f"关联的 Site 不存在: {row['site']}")


def add_MAG():
    df = pd.read_excel(f'db/static/db/MAG.xlsx')

    for index, row in df.iterrows():
        try:
            # 获取关联的 Site 对象
            sample = Sample.objects.get(id=row['sample_id'])

            # 创建 Sample 对象
            mm = MAG(
                id=row['id'],
                sample=sample,
                Kingdom=row.get('Kingdom'),
                Phylum=row.get('Phylum'),
                Class=row.get('Class'),
                Order=row.get('Order'),
                Family=row.get('Family'),
                Genus=row.get('Genus'),
                Species=row.get('Species'),
                taxonomy=row.get('taxonomy'),
                relative_abundance=row.get('relative_abundance'),
                length=row.get('length'),
                RPKM=row.get('RPKM'),
                completeness=row.get('completeness'),
                contamination=row.get('contamination'),
                GC=row.get('GC'),
                N50 = row.get('N50'),
                old=row.get('name'),
            )

            mm.save()
            print(f"MAG {row['id']} 已成功保存！")
        except Site.DoesNotExist:
            print(f"关联的 Site 不存在: {row['site']}")



def export_sample_mag_taxonomy():

    samples = Sample.objects.all()  # 修改为实际的 Sample ID
    for sample in samples:
        # 获取该 Sample 关联的 MAG 的 taxonomy 信息
        mags = MAG.objects.filter(sample=sample).values('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

        # 将数据转为 pandas DataFrame
        df = pd.DataFrame(list(mags))
        df = df.replace({'': '-'})
        # 定义 Excel 文件的名称
        filename = r'G:\lab\db\static\db\sample_tax/'+sample.id+'.xlsx'


        df.to_excel(filename, index=False,header=False)

def export_site_mag_taxonomy():

    sites = Site.objects.all()  # 修改为实际的 Sample ID
    for site in sites:
        # 获取该 Sample 关联的 MAG 的 taxonomy 信息
        mags = MAG.objects.filter(sample__site=site).values('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species')

        # 将数据转为 pandas DataFrame
        df = pd.DataFrame(list(mags))

        # 定义 Excel 文件的名称
        filename = r'G:\lab\db\static\db\site_tax/'+site.id+'.xlsx'


        df.to_excel(filename, index=False,header=False)

def export_pathway_step_radio():
    samples = Sample.objects.all()
    for sample in samples:
        s = {}
        m = []
        mags = MAG.objects.filter(sample=sample)
        for mag in mags:
            m.append(mag.id)
            paths=Pathway.objects.filter(MAG=mag)
            p = {}
            for path in paths:
                p[path.module_ID] = round(path.steps_percentage,2)
            s[mag.id] = p
        # transposed_list = [[s[j][i] for j in range(len(s))] for i in range(len(s[0]))]

        result = [
            [index, value, key]
            for index, (rw_id, pathways) in enumerate(s.items())
            for key, value in pathways.items()
        ]
        data = {'mags': m,'data':result}
        output_file = r'G:\lab\db\static\db\river\\'+sample.id+'.json'
        with open(output_file, "w") as f:
            json.dump(data, f, indent=4)  # indent=4 用于更好地格式化输出

        print(f"列表已保存为 JSON 文件: {output_file}")

def extract_edges(node, parent=None):
    edges = []
    if "name" in node:
        if parent:  # 如果有父节点
            edge = {
                "source": parent,
                "target": node["name"],
                "value": node.get("value", None)  # 子节点的值
            }
            edges.append(edge)
        # 递归处理子节点
        if "children" in node:
            for child in node["children"]:
                edges.extend(extract_edges(child, node["name"]))
    return edges

def rename_id():
    folder_path = r"G:\lab\db\static\db\annotation"
    for file in os.listdir(folder_path):
        if file.endswith(".annotations"):
            # 提取文件名前半部分（去掉 .fasta 后缀和其他后缀部分）
            base_name = file.split(".")[0] + '.' + file.split(".")[1]  # 获取第一个 `.` 之前的部分

            a = MAG.objects.get(old=base_name)
            new_name = a.id + '.annotations'
            old_path = os.path.join(folder_path, file)
            new_path = os.path.join(folder_path, new_name)

            os.rename(old_path, new_path)
            print(f"重命名: {file} → {new_name}")


    print("文件重命名完成！")

def ModifyGeneInFasta():
    # 定义需要处理的文件夹路径
    folder_path = r"G:\lab\db\static\db\fasta"

    # 遍历文件夹中的所有 .fa 文件
    for filename in os.listdir(folder_path):

        file_path = os.path.join(folder_path, filename)
        new_lines = []

        # 打开并读取 .fa 文件
        with open(file_path, "r") as file:
            for line in file:

                if line.startswith(">"):
                    # 在序列名称前加上文件名（不包括扩展名）
                    id = Gene.objects.get(name=line[1:].strip()).id

                    new_line = f">{id}{line[-1]}"
                    new_lines.append(new_line)
                else:
                    new_lines.append(line)
        print(file_path)

        # 将修改后的内容写回文件，覆盖原文件
        with open(file_path, "w") as file:
            file.writelines(new_lines)
    print("所有 .fa 文件的标题已成功更新。")


def ModifyGeneInAnnotation():
    # 定义需要处理的文件夹路径
    folder_path = r"G:\lab\db\static\db\annotation"
    # 遍历文件夹中的所有 .fa 文件
    for filename in os.listdir(folder_path):
        file_path = os.path.join(folder_path, filename)
        new_lines = []

        with open(file_path, 'r') as f:
            for line in f:
                # 拆分行数据
                if line.startswith('#'):
                    new_lines.append(line)
                    continue
                columns = line.strip().split('\t')

                name = columns[0]
                columns[0] = Gene.objects.get(name=name).id
                new_line = "\t".join(columns) + "\n"
                new_lines.append(new_line)

        # 将修改后的内容写回文件，覆盖原文件
        with open(file_path, "w") as file:
            file.writelines(new_lines)
        print(name)
    print("所有 .annotation 文件的标题已成功更新。")


def export_sankey():

    sites = Site.objects.all()
    for site in sites:

        nodes = set()
        links = []
        samples = Sample.objects.filter(site=site)
        for sample in samples:
            nodes.add(sample.id)

            mags = MAG.objects.filter(sample=sample)

            ma = MAG.objects.filter(Q(sample=sample) & Q(Kingdom='Archaea'))
            mb = MAG.objects.filter(Q(sample=sample) & Q(Kingdom='Bacteria'))
            if len(ma) != 0:
                links.append({'source': sample.id, 'target': 'k_Archaea', 'value': len(ma)})
            links.append({'source': sample.id, 'target': 'k_Bacteria', 'value': len(mb)})
            for mag in mags:
                nodes.add('k_'+mag.Kingdom)
                nodes.add('p_'+mag.Phylum)
                nodes.add('c_'+mag.Class)
                nodes.add('o_'+mag.Order)
                nodes.add('f_'+mag.Family)
                nodes.add('g_'+mag.Genus)
                nodes.add('s_'+mag.Species)

        input_path = r"G:\db_data_process\sunburst_data\site_json/" + site.id + ".json"
        with open(input_path, "r", encoding="utf-8") as f:
            data = json.load(f)  # 加载 JSON 数据
        edges = extract_edges(data)
        edges.pop(0)
        links.extend(edges)
        node_list = []
        for node in nodes:
            node_list.append({"name": node})
        data = {'nodes': node_list, 'links': links}
        print(len(node_list))
        output_file = r'G:\lab\db\static\db\sankey\\'+site.id+'.json'
        with open(output_file, "w") as f:
            json.dump(data, f, indent=4)  # indent=4 用于更好地格式化输出
def notknow():
    geneDic = {'PRK': 'K00855', 'rbcL': 'K01601', 'rbcS': 'K01602', 'PGK': 'K00927', 'GAPA': 'K05298', 'gap2': 'K00150',
               'GAPDH': 'K00134', 'TPI': 'K01803', 'ALDO': 'K01623', 'FBA': 'K01624', 'FBP': 'K03841', 'glpX': 'K02446',
               'glpX-SEBP': 'K11532', 'fbp-SEBP': 'K01086', 'E2.2.1.1': 'K00615', 'E3.1.3.37': 'K01100',
               'rpiA': 'K01807', 'rpiB': 'K01808', 'rpe': 'K01783',
               'porA': 'K00169', 'porB': 'K00170', 'porD': 'K00171', 'porC': 'K00172', 'por': 'K03737', 'pps': 'K01007',
               'ppdK': 'K01006', 'ppc': 'K01595', 'pycA': 'K01959', 'pycB': 'K01960', 'PC': 'K01958', 'mdh': 'K00024',
               'E4.2.1.2A': 'K01676', 'E4.2.1.2B': 'K01679',
               'E4.2.1.2AA': 'K01677', 'E4.2.1.2AB': 'K01678', 'sdhA': 'K00239', 'sdhB': 'K00240', 'sdhC': 'K00241',
               'sdhD': 'K00242', 'frdA': 'K18556', 'frdB': 'K18557', 'frdC': 'K18558', 'frdD': 'K18559',
               'frdE': 'K18560', 'sucD': 'K01902', 'sucC': 'K01903', 'korA': 'K00174',
               'korB': 'K00175', 'korC': 'K00177', 'korD': 'K00176', 'IDH1': 'K00031', 'ACO': 'K01681',
               'acnA': 'K27802', 'acnB': 'K01682', 'aclA': 'K15230', 'aclB': 'K15231', 'ccsA': 'K15232',
               'ccsB': 'K15233', 'ccl': 'K15234', 'K15038': 'K15038', 'K15017': 'K15017', 'K14465': 'K14465',
               '4hbl': 'K14467',
               'K18861': 'K18861', 'K25774': 'K25774', 'abfD': 'K14534', 'K15016': 'K15016', 'ACAT': 'K00626',
               'K01964': 'K01964', 'K15037': 'K15037', 'K15036': 'K15036', 'K15039': 'K15039', 'K15018': 'K15018',
               'K15019': 'K15019', 'K15020': 'K15020', 'MCEE': 'K05606',
               'E5.4.99.2A': 'K01848', 'E5.4.99.2B': 'K01849', 'K14466': 'K14466', 'accB': 'K02160', 'accC': 'K01961',
               'accA': 'K01962', 'accD': 'K01963', 'mcr': 'K14468', 'K14469': 'K14469', 'K15052': 'K15052',
               'MUT': 'K01847', 'smtA1': 'K14471',
               'smtB': 'K14472', 'mcl': 'K08691', 'mch': 'K14449', 'mct': 'K14470', 'meh': 'K09709', 'cooS': 'K00198',
               'fdhA': 'K05299', 'fdhB': 'K15022', 'fdhF': 'K22015', 'hydA2': 'K25123', 'hycB': 'K25124',
               'fhs': 'K01938', 'folD': 'K01491', 'fchA': 'K01500', 'metF': 'K00297', 'metV': 'K25007',
               'rnfC2': 'K25008', 'acsE': 'K15023', 'acsB': 'K14138', 'cdhE': 'K00197', 'cdhD': 'K00194'
               }
    module_map = {
        'M00165': [
            {'name': 'Step 1', 'data': 1, 'genes': ['PRK']},
            {'name': 'Step 2', 'data': 2, 'genes': ['rbcL', 'rbcS']},
            {'name': 'Step 3', 'data': 1, 'genes': ['PGK']},
            {'name': 'Step 4', 'data': 3, 'genes': ['GAPA', 'gap2', 'GAPDH']},
            {'name': 'Step 5', 'data': 1, 'genes': ['TPI']},
            {'name': 'Step 6', 'data': 2, 'genes': ['ALDO', 'FBA']},
            {'name': 'Step 7', 'data': 4, 'genes': ['FBP', 'glpX', 'glpX-SEBP', 'fbp-SEBP']},
            {'name': 'Step 8', 'data': 1, 'genes': ['E2.2.1.1']},
            {'name': 'Step 9', 'data': 3, 'genes': ['E3.1.3.37', 'glpX-SEBP', 'fbp-SEBP']},
            {'name': 'Step 10', 'data': 2, 'genes': ['rpiA', 'rpiB']},
            {'name': 'Step 11', 'data': 1, 'genes': ['rpe']},
        ],
        'M00173': [
            {'name': 'Step 1', 'data': 5, 'genes': ['porA', 'porB', 'porD', 'porC', 'por']},
            {'name': 'Step 2', 'data': 6, 'genes': ['pps', 'ppdK', 'ppc', 'pycA', 'pycB', 'PC']},
            {'name': 'Step 3', 'data': 1, 'genes': ['mdh']},
            {'name': 'Step 4', 'data': 4, 'genes': ['E4.2.1.2A', 'E4.2.1.2B', 'E4.2.1.2AA', 'E4.2.1.2AB']},
            {'name': 'Step 5', 'data': 9,
             'genes': ['sdhA', 'sdhB', 'sdhC', 'sdhD', 'frdA', 'frdB', 'frdC', 'frdD', 'frdE']},
            {'name': 'Step 6', 'data': 2, 'genes': ['sucD', 'sucC']},
            {'name': 'Step 7', 'data': 4, 'genes': ['korA', 'korB', 'korC', 'korD']},
            {'name': 'Step 8', 'data': 1, 'genes': ['IDH1']},
            {'name': 'Step 9', 'data': 3, 'genes': ['ACO', 'acnA', 'acnB']},
            {'name': 'Step 10', 'data': 5, 'genes': ['aclA', 'aclB', 'ccsA', 'ccsB', 'ccl']},
        ],
        'M00374': [
            {'name': 'Step 1', 'data': 4, 'genes': ['porA', 'porB', 'porD', 'porC']},
            {'name': 'Step 2', 'data': 1, 'genes': ['pps']},
            {'name': 'Step 3', 'data': 1, 'genes': ['ppc']},
            {'name': 'Step 4', 'data': 1, 'genes': ['mdh']},
            {'name': 'Step 5', 'data': 2, 'genes': ['E4.2.1.2AA', 'E4.2.1.2AB']},
            {'name': 'Step 6', 'data': 4, 'genes': ['sdhA', 'sdhB', 'sdhC', 'sdhD']},
            {'name': 'Step 7', 'data': 2, 'genes': ['sucD', 'sucC']},
            {'name': 'Step 8', 'data': 2, 'genes': ['K15038', 'K15017']},
            {'name': 'Step 9', 'data': 1, 'genes': ['K14465']},
            {'name': 'Step 10', 'data': 3, 'genes': ['4hbl', 'K18861', 'K25774']},
            {'name': 'Step 11', 'data': 1, 'genes': ['abfD']},
            {'name': 'Step 12', 'data': 1, 'genes': ['K15016']},
            {'name': 'Step 13', 'data': 1, 'genes': ['ACAT']},
        ],
        'M00375': [
            {'name': 'Step 1', 'data': 3, 'genes': ['K01964', 'K15037', 'K15036']},
            {'name': 'Step 2', 'data': 1, 'genes': ['K15017']},
            {'name': 'Step 3', 'data': 1, 'genes': ['K15039']},
            {'name': 'Step 4', 'data': 1, 'genes': ['K15018']},
            {'name': 'Step 5', 'data': 1, 'genes': ['K15019']},
            {'name': 'Step 6', 'data': 1, 'genes': ['K15020']},
            {'name': 'Step 7', 'data': 1, 'genes': ['MCEE']},
            {'name': 'Step 8', 'data': 2, 'genes': ['E5.4.99.2A', 'E5.4.99.2B']},
            {'name': 'Step 9', 'data': 2, 'genes': ['K15038', 'K15017']},
            {'name': 'Step 10', 'data': 1, 'genes': ['K14465']},
            {'name': 'Step 11', 'data': 3, 'genes': ['K14466', 'K18861', 'K25774']},
            {'name': 'Step 12', 'data': 1, 'genes': ['abfD']},
            {'name': 'Step 13', 'data': 1, 'genes': ['K15016']},
            {'name': 'Step 14', 'data': 1, 'genes': ['ACAT']},
        ],
        'M00376': [
            {'name': 'Step 1', 'data': 4, 'genes': ['accB', 'accC', 'accA', 'accD']},
            {'name': 'Step 2', 'data': 1, 'genes': ['mcr']},
            {'name': 'Step 3', 'data': 1, 'genes': ['K14469']},
            {'name': 'Step 4', 'data': 1, 'genes': ['K15052']},
            {'name': 'Step 5', 'data': 1, 'genes': ['MCEE']},
            {'name': 'Step 6', 'data': 3, 'genes': ['MUT', 'E5.4.99.2A', 'E5.4.99.2B']},
            {'name': 'Step 7', 'data': 2, 'genes': ['smtA1', 'smtB']},
            {'name': 'Step 8', 'data': 3, 'genes': ['sdhA', 'sdhB', 'sdhC']},
            {'name': 'Step 9', 'data': 1, 'genes': ['E4.2.1.2B']},
            {'name': 'Step 10', 'data': 1, 'genes': ['mcl']},
            {'name': 'Step 11', 'data': 1, 'genes': ['mch']},
            {'name': 'Step 12', 'data': 1, 'genes': ['mct']},
            {'name': 'Step 13', 'data': 1, 'genes': ['meh']},
        ],
        'M00377': [
            {'name': 'Step 1', 'data': 1, 'genes': ['cooS']},
            {'name': 'Step 2', 'data': 5, 'genes': ['fdhA', 'fdhB', 'fdhF', 'hydA2', 'hycB']},
            {'name': 'Step 3', 'data': 1, 'genes': ['fhs']},
            {'name': 'Step 4', 'data': 2, 'genes': ['folD', 'fchA']},
            {'name': 'Step 5', 'data': 3, 'genes': ['metF', 'metV', 'rnfC2']},
            {'name': 'Step 6', 'data': 1, 'genes': ['acsE']},
            {'name': 'Step 7', 'data': 3, 'genes': ['acsB', 'cdhE', 'cdhD']},
        ]
    }
    for module in module_map:
        genes = {}
        for step in module_map[module]:
            for i in step['genes']:
                genes[geneDic[i]] = i
        print(genes)

def getAutoCompleteData():
    mags = MAG.objects.all()
    l = set()
    for mag in mags:
        l.add(mag.Kingdom)
        l.add(mag.Phylum)
        l.add(mag.Class)
        l.add(mag.Order)
        l.add(mag.Family)
        l.add(mag.Genus)
        l.add(mag.Species)
    content = "\n".join(map(str, l))

    with open("getAutoCompleteData.txt", "w", encoding="utf-8") as f:
        f.write(content)
if __name__ == "__main__":
    import pandas as pd


    # 1. 读取 Excel 文件
    df = pd.read_excel(r"G:\lab\db\static\db\pathway\M00377.xlsx")

    # 2. 获取所有 MAG.id 和 name 映射成字典
    # 假设 Excel 第一列为 MAG.id，对应要替换为 MAG.name
    id_list = df[df.columns[0]].tolist()  # 获取第一列所有 ID

    # 批量从数据库中获取对应 name
    mag_map = {str(mag.old): mag.id for mag in MAG.objects.filter(old__in=id_list)}

    # 3. 替换 Excel 第一列
    df[df.columns[0]] = df[df.columns[0]].astype(str).map(mag_map).fillna(df[df.columns[0]])

    # 4. 保存新 Excel 文件
    df.to_excel("M00377.xlsx", index=False)







