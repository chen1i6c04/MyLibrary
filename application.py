import base64
import requests
import subprocess
from Bio.Application import AbstractCommandline, _Option, _Switch


class MlstCommandline(AbstractCommandline):
    def __init__(self, cmd='mlst.py', **kwargs):
        self.parameters = [
            _Option(['-i', 'infile'], '', equate=False, is_required=True),
            _Option(['-o', 'outdir'], '', equate=False, is_required=True),
            _Option(['-p', 'database'], '', equate=False, is_required=True),
            _Option(['-s', 'species'], '', equate=False, is_required=True),
            _Option(['-t', 'tmp'], '', equate=False),
            _Option(['-mp', 'method_path'], '', equate=False),
            _Switch(['-x', 'extented_output'], ''),
            _Switch(['-q', 'quiet'], ''),
            _Switch(['-matrix', 'matrix'], '')
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class PlasmidFinderCommandline(AbstractCommandline):
    def __init__(self, cmd='plasmidfinder.py', **kwargs):
        self.parameters = [
            _Option(['-i', 'infile'], '', equate=False, is_required=True),
            _Option(['-o', 'outdir'], '', equate=False, is_required=True),
            _Option(['-p', 'database'], '', equate=False, is_required=True),
            _Option(['-tmp', 'tmp'], '', equate=False),
            _Option(['-l', 'mincov'], '', equate=False),
            _Option(['-t', 'threshold'], '', equate=False),
            _Switch(['-x', 'extented_output'], ''),
            _Switch(['-q', 'quiet'], ''),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class VirulencefinderCommandline(AbstractCommandline):
    def __init__(self, cmd='virulencefinder.py', **kwargs):
        self.parameters = [
            _Option(['-i', 'infile'], '', equate=False, is_required=True),
            _Option(['-o', 'outdir'], '', equate=False, is_required=True),
            _Option(['-p', 'db_path'], '', equate=False, is_required=True),
            _Option(['-d', 'database'], '', equate=False, is_required=True),
            _Option(['-tmp', 'tmp'], '', equate=False),
            _Option(['-mp', 'method_path'], '', equate=False),
            _Option(['--mincov', 'mincov'], '', equate=False),
            _Option(['--threshold', 'threshold'], '', equate=False),
            _Switch(['-x', 'extented_output'], ''),
            _Switch(['-q', 'quiet'], ''),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class ResfinderCommandline(AbstractCommandline):
    def __init__(self, cmd='run_resfinder.py', **kwargs):
        self.parameters = [
            _Option(['-ifa', 'infasta'], '', equate=False, is_required=False),
            _Option(['-ifq', 'infastq'], '', equate=False, is_required=False),
            _Option(['-o', 'outdir'], '', equate=False, is_required=True),
            _Option(['--db_path_point', 'db_point'], '', equate=False, is_required=True),
            _Option(['--db_path_res', 'db_res'], '', equate=False, is_required=True),
            _Option(['--species', 'species'], '', equate=False, is_required=False),
            _Option(['--min_cov', 'min_cov'], '', equate=False, is_required=False),
            _Option(['--threshold', 'threshold'], '', equate=False, is_required=False),
            _Option(['--min_cov_point', 'min_cov_point'], '', equate=False, is_required=False),
            _Option(['--threshold_point', 'threshold_point'], '', equate=False, is_required=False),
            _Switch(['-u', 'unknown'], ''),
            _Switch(['-c', 'point'], ''),
            _Switch(['-acq', 'acquired'], ''),
            _Switch(['--nanopore', 'nanopore'], ''),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class SerotypefinderCommandline(AbstractCommandline):
    def __init__(self, cmd='serotypefinder.py', **kwargs):
        self.parameters = [
            _Option(['-i', 'infile'], '', equate=False, is_required=True),
            _Option(['-o', 'outdir'], '', equate=False, is_required=True),
            _Option(['-p', 'db_path'], '', equate=False, is_required=True),
            _Option(['-d', 'database'], '', equate=False, is_required=False),
            _Option(['-tmp', 'tmp'], '', equate=False),
            _Option(['-mp', 'method_path'], '', equate=False),
            _Switch(['-x', 'extented_output'], ''),
            _Switch(['-q', 'quiet'], ''),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class AmrfinderCommandline(AbstractCommandline):
    def __init__(self, cmd="amrfinder", **kwargs):
        self.parameters = [
            _Option(['-n', 'nuc_fasta'], '', equate=False, is_required=True),
            _Option(['-o', 'output_file'], '', equate=False, is_required=True),
            _Option(['-d', 'database'], '', equate=False, is_required=True),
            _Option(['-O', '--organism', 'organism'], '', equate=False, is_required=False),
            _Option(['--nucleotide_output', 'nucleotide_output'], '', equate=False, is_required=False),
            _Option(['--threads', 'threads'], '', equate=False, is_required=False),
            _Switch(['--mutation_all', 'mutation_all'], ''),
            _Switch(['--plus', 'plus'], ''),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


def ribosomal_mlst(upload_file, output):
    uri = 'http://rest.pubmlst.org/db/pubmlst_rmlst_seqdef_kiosk/schemes/1/sequence'
    with open(upload_file, 'r') as x:
        fasta = x.read()
    payload = '{"base64":true,"details":true,"sequence":"' + base64.b64encode(fasta.encode()).decode() + '"}'
    response = requests.post(uri, data=payload)
    if response.status_code == requests.codes.ok:
        data = response.json()
        try:
            with open(output, 'w') as f:
                json.dump(data['taxon_prediction'], f)
        except ValueError:
            pass


class ProdigalCommandline(AbstractCommandline):
    """
    :param kwargs:
        a:  Write protein translations to the selected file.
        c:  Closed ends.  Do not allow genes to run off edges.
        d:  Write nucleotide sequences of genes to the selected file.
        f:  Select output format (gbk, gff, or sco).  Default is gbk.
        g:  Specify a translation table to use (default 11).
        h:  Print help menu and exit.
        i:  Specify FASTA/Genbank input file (default reads from stdin).
        m:  Treat runs of N as masked sequence; don't build genes across them.
        n:  Bypass Shine-Dalgarno trainer and force a full motif scan.
        o:  Specify output file (default writes to stdout).
        p:  Select procedure (single or meta).  Default is single.
        q:  Run quietly (suppress normal stderr output).
        s:  Write all potential genes (with scores) to the selected file.
        t:  Write a training file (if none exists); otherwise, read and use
            the specified training file.
        v:  Print version number and exit.
    :return:
        str: prodigal command line
    """
    def __init__(self, cmd="prodigal", **kwargs):
        self.parameters = [
            _Option(['-a', 'a'], '', equate=False, is_required=False),
            _Option(['-c', 'c'], '', equate=False, is_required=False),
            _Option(['-d', 'd'], '', equate=False, is_required=False),
            _Option(['-f', 'f'], '', equate=False, is_required=False),
            _Option(['-g', 'g'], '', equate=False, is_required=False),
            _Option(['-i', 'i'], '', equate=False, is_required=True),
            _Option(['-m', 'm'], '', equate=False, is_required=False),
            _Option(['-n', 'n'], '', equate=False, is_required=False),
            _Option(['-o', 'o'], '', equate=False, is_required=False),
            _Option(['-p', 'p'], '', equate=False, is_required=False),
            _Option(['-q', 'q'], '', equate=False, is_required=False),
            _Option(['-s', 's'], '', equate=False, is_required=False),
            _Option(['-t', 't'], '', equate=False, is_required=False),
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)


class FastANICommandline(AbstractCommandline):
    def __init__(self, cmd='fastANI', **kwargs):
        self.parameters = [
            _Option(['-r', 'r'], '', equate=False),
            _Option(['--rl', 'rl'], '', equate=False),
            _Option(['-q', 'q'], '', equate=False),
            _Option(['--ql', 'ql'], '', equate=False),
            _Option(['-o', 'output'], '', equate=False),
            _Option(['-k', 'kmer'], '', equate=False),
            _Option(['-t', 'threads'], '', equate=False)
        ]
        AbstractCommandline.__init__(self, cmd, **kwargs)
        

def assembly_stats(input_file):
    output = subprocess.getoutput(f'assembly-stats -t {input_file}')
    head, values = list(map(lambda x: x.split(), output.splitlines()))
    return dict(zip(head, values))
