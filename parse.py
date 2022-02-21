import json
import pandas as pd


def parse_amrfinder_result(file):
    fields = {
        'Contig id': 'contig_id',
        'Gene symbol': 'gene_symbol',
        'Element subtype': 'element_subtype',
        'Method': 'method',
        '% Coverage of reference sequence': 'coverage_of_reference_sequence',
        '% Identity to reference sequence': 'identity_of_reference_sequence'}
    records = pd.read_csv(file, sep='\t', usecols=fields.keys()).rename(columns=fields).to_dict(orient='records')
    return records


def parse_mlst_result(file):
    record = dict()
    with open(file) as f:
        data = json.load(f)
    sequence_type = data['mlst']['results']['sequence_type']
    sequence_type = sequence_type if sequence_type.isnumeric() else ''
    record['ST'] = sequence_type
    allele_profile = data['mlst']['results']['allele_profile']
    record.update({allele_name: (profile['allele'] if profile['allele'].isnumeric() else '')
                   for allele_name, profile in allele_profile.items()})
    return record


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