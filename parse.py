import re
import os
import json
import pandas as pd


def parse_amrfinder_result(file):
    fields = {
        'Contig id': 'contig_id',
        'Gene symbol': 'gene_symbol',
        'Element subtype': 'element_subtype',
        'Method': 'method',
        '% Coverage of reference sequence': 'coverage_of_reference_sequence',
        '% Identity to reference sequence': 'identity_of_reference_sequence',
    }
    df = pd.read_csv(file, sep='\t', usecols=fields.keys()).rename(columns=fields)
    return {subtype: group.to_dict(orient='records') for subtype, group in df.groupby('element_subtype')}


def parse_mlst_result(file):
    with open(file) as handle:
        data = json.load(handle)[0]
    return data


def parse_plasmidfinder_result(file):
    df = pd.read_csv(file, sep='\t', usecols=['Plasmid'])
    return df['Plasmid'].drop_duplicates().to_list()


def parse_sistr_result(file):
    with open(file) as f:
        data = json.load(f)[0]
    summary = dict()
    summary['subspecies'] = data['cgmlst_subspecies']
    summary['serovar'] = data['serovar']
    summary['serogroup'] = data['serogroup']
    return summary


def parse_lissero_result(file):
    with open(file) as f:
        next(f)
        return f.read().split('\t')[1]


def parse_virulencefinder_result(file):
    with open(file) as f:
        data = json.load(f)
    hits = data['virulencefinder']['results']['Listeria']['listeria']
    return sorted(set(i['virulence_gene'] for i in hits.values()))


def parse_pointfinder_result(file):
    df = pd.read_csv(file, sep='\t', usecols=['Mutation'])
    return df['Mutation'].sort_values().str.replace(' p.', '_').str.replace(' r.', '_').to_list()


def parse_resfinder_result(file):
    df = pd.read_csv(file, sep='\t', usecols=['Resistance gene'])
    return df['Resistance gene'].drop_duplicates().sort_values().str.cat(sep=', ')


def parse_busco_summary(file):
    lieage = ""
    complete = ""
    singel = ""
    duplicated = ""
    fragmented = ""
    missing = ""
    with open(file) as handle:
        for line in handle.readlines():
            result = re.search(': ([a-z]+_odb[0-9]+) ', line)
            if result:
                lieage = result.groups()[0]   
            result = re.findall('[A-Z]:([0-9]+.[0-9])%', line)
            if result:
                complete, singel, duplicated, fragmented, missing = [float(i) for i in result]
    return lieage, complete, singel, duplicated, fragmented, missing

def parse_busco_summary2(file):
    with open(file) as handle:
        data = json.load(handle)
    total_markers = int(data['dataset_total_buscos'])
    lineage = os.path.basename(data['dataset'])
    complete = round(data['C']/total_markers*100, 1)
    single = round(data['S']/total_markers*100, 1)
    duplicated = round(data['D']/total_markers*100, 1)
    fragmented = round(data['F']/total_markers*100, 1)
    missing = round(data['M']/total_markers*100, 1)
    return lineage, complete, single, duplicated, fragmented, missing


def parse_fastani_result(file):
    with open(file) as handle:
        return float(next(handle).strip().split()[2])
        