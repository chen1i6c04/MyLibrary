import os
import shutil
import subprocess
from tempfile import TemporaryDirectory
from application import PlasmidFinderCommandline, ResfinderCommandline, VirulencefinderCommandline,\
    SerotypefinderCommandline, AmrfinderCommandline


def run_mlst(infile, output_json, output_fasta, scheme=None, threads=2):
    cmd = ['mlst', '--json', output_json, '--novel', output_fasta, '--nopath', infile, '--threads', str(threads)]
    if scheme:
        cmd += ['--scheme', scheme]
    p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    if p.returncode != 0:
        raise subprocess.CalledProcessError(p.returncode, cmd, output=p.stdout, stderr=p.stderr)


def run_plasmidfinder(infile, outdir, database, identity=0.6, mincov=0.9):
    os.makedirs(outdir, exist_ok=True)
    with TemporaryDirectory(dir='/dev/shm/') as tmp:
        cline = PlasmidFinderCommandline(
            infile=infile, outdir=outdir, database=database, tmp=tmp, extented_output=True, mincov=mincov, threshold=identity
        )
        cline()


def run_resfinder(infile, outdir, db_res, db_point, **kwargs):
    """
    :param kwargs:
    species
    unknown
    point
    acquired
    min_cov
    threshold
    """
    os.makedirs(outdir, exist_ok=True)
    cline = ResfinderCommandline(
        infasta=infile, outdir=outdir, db_res=db_res, db_point=db_point, **kwargs
    )
    cline()
    for dirname in ('pointfinder_blast', 'resfinder_blast'):
        if os.path.exists(os.path.join(outdir, dirname)):
            shutil.rmtree(os.path.join(outdir, dirname))


def run_amrfinder(infile, outfile, database, **kwargs):
    """
    :param kwargs:
    organism
    nucleotide_output
    threads
    mutation_all
    """
    cline = AmrfinderCommandline(
        cmd='/home/chen1i6c04/miniconda3/envs/amrfinder/bin/amrfinder',
        nuc_fasta=infile, output_file=outfile, database=database, **kwargs
    )
    cline()


def run_virulencefinder(infile, outdir, db_path, database):
    """
    database: virulence_ecoli, virulence_ent, listeria, s.aureus_exoenzyme, s.aureus_hostimm, s.aureus_toxin, stx
    :return:
    """
    os.makedirs(outdir, exist_ok=True)
    with TemporaryDirectory(dir='/tmp/') as tmpdir:
        cline = VirulencefinderCommandline(
            infile=infile, outdir=outdir, db_path=db_path, database=database, tmp=tmpdir, extented_output=True
        )
        cline()