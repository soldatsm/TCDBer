import pandas as pd
import subprocess  # for bash scripts
import argparse
import shutil  # check installed programs
import sys  # check installed modules
import requests as re
import json
from tqdm import tqdm
from logo import logo

parser = argparse.ArgumentParser()

parser.add_argument('-check_modules', help='Check required packages and programmes')
parser.add_argument('-make_query', nargs='+', help='If you have table with names and '
                                                   'AA sequences you can make fasta file for blast.'
                                                   'In addition you can specifies table delimiter. '
                                                   'This command accepts 1 - path to table, '
                                                   '2 - index of protein ids column, '
                                                   '3 - index of protein seq column, '
                                                   '4 - FASTA file name,'
                                                   '5 - delimiter (def=comma)')
parser.add_argument('-download_TCDB', help='Download TCDB database')
parser.add_argument('-make_blast_db', help='Make blastp db from TCDB data')
parser.add_argument('-blastp', nargs='+', help='Do blastp enter query name and name of output tsv '
                                               'table separated by space')
parser.add_argument('-make_table', nargs='+',
                    help='make table form blast output, input - path to blast output and output file.tsv')
parser.add_argument('-evidences', help='Add evidence information to the table')

args = parser.parse_args()


def module_checker():
    require_mod_lst = ['pandas']
    require_prog_lst = ['makeblastdb', 'blastp', 'wget']
    cheker = True
    print('Checking \n')
    for prog in require_prog_lst:
        if shutil.which(prog) is not None:
            print(f'{prog} -- installed')
        else:
            cheker = False
            print(f'{prog} -- is not installed')

    for mod in require_mod_lst:
        if mod in sys.modules:
            print(f'{mod} -- installed')
        else:
            cheker = False
            print(f'{mod} -- is not installed')
    if cheker is False:
        print('\n *Please install missed pakages or programs* \n')
        sys.exit()
    else:
        print('\n')
        print('-------------------------------------------')
        print(r'All required modules\programs are installed')
        print('-------------------------------------------\n')


def query(table: str,
          name_idx: int,
          seq_idx: int,
          output_name: str, delim=','):

    tab = pd.read_csv(table, delimiter=delim)
    with open(output_name, 'w+') as write_file:
        for idx, row in tab.iterrows():
            write_file.write('>' +
                             row.iloc[name_idx] +
                             '\n' +
                             row.iloc[seq_idx] +
                             '\n')
    print('Done')


def download_db():
    with open(f'./db_loader.sh', 'w') as write_file:
        write_file.write(f"""\
        #!/bin/sh
        wget -O tcdb_db_"$(date +%Y-%m-%d)".fasta http://www.tcdb.org/public/tcdb
        
        echo "TCDB data download completed"
        """)
    subprocess.run(['chmod', '+x', './db_loader.sh'])
    subprocess.run(['./db_loader.sh'], shell=True)


def make_db():
    with open(f'./blastp_db_maker.sh', 'w') as write_file:
        write_file.write(f"""\
        #!/bin/sh
        makeblastdb -in ./tcdb_db* -dbtype prot -title tcdb_db_"$(date)"
        """)
        subprocess.run(['chmod', '+x', './blastp_db_maker.sh'])


def make_blastp(query: str,
                blast_results_name: str,
                db_name='tcdb_db_*.fasta',
                threads=1,
                script_name='blastp_script.sh'):
    """
    :param query:
    :param db_name:
    :param threads:
    :param script_name:
    :param blast_results_name:
    :return: tsv file, blast.sh

    qaccver	saccver	stitle	pident	length	mismatch	gapopen	evalue	qcovs	bitscore
    tabdelimited

    This function make blastp script and execute it
    """

    with open(f'./{script_name}', 'w') as write_file:
        write_file.write(f"""\
        #!/bin/sh
        blastp -query {query} -db {db_name} -num_threads {threads} -outfmt "6 qaccver stitle pident length mismatch gapopen evalue qcovs bitscore" -out {blast_results_name}
        """)

    subprocess.run(['chmod', '+x', f'./{script_name}'])  # since script is unexecutable after creating
    subprocess.run([f'./{script_name}'], shell=True)


def remake_table(blast_out_tab: str, output_name: str) -> None:
    data = pd.read_csv(blast_out_tab, delimiter='\t',
                       header=None)
    qaccver = []
    descipion = []
    TCDB_id = []
    UniprotID = []
    pident = []
    length = []
    mismatch = []
    gapopen = []
    evalue = []
    qcovs = []
    bitscore = []

    for idx, row in data.iterrows():
        descipion_s = row.iloc[1].split('OS')[0].split(' ')[1::]
        descipion_s = ' '.join(descipion_s)
        name_s = row.iloc[1].split(' ')[0]
        TCDB_id.append(name_s.split('|')[3])
        UniprotID.append(name_s.split('|')[2])
        qaccver.append(row.iloc[0])
        descipion.append(descipion_s)
        pident.append(row.iloc[2])
        length.append(row.iloc[3])
        mismatch.append(row.iloc[4])
        gapopen.append(row.iloc[5])
        evalue.append(row.iloc[6])
        qcovs.append(row.iloc[7])
        bitscore.append(row.iloc[8])

    output_table = pd.DataFrame({
        'qaccver': qaccver,
        'UiprotID': UniprotID,
        'TCDB_id': TCDB_id,
        'description': descipion,
        'pident': pident,
        'length': length,
        'mismatch': mismatch,
        'gapopen': gapopen,
        'evalue': evalue,
        'qcovs': qcovs,
        'bitscore': bitscore})

    output_table.to_csv(output_name, sep='\t', index=False)


def evidanceer(table, col_idx=1) -> None:
    col_idx: int
    tab = pd.read_csv(table, delimiter='\t')

    dedup_table = tab#.drop_duplicates(keep='first', subset='qaccver')

    json_lst = []
    urls = []

    for idx, row in tqdm(dedup_table.iterrows()):
        url = f'https://rest.uniprot.org/uniprotkb/search?query={row.iloc[col_idx]}'
        urls.append(url)
        if re.get(url):
            handler = re.get(url)
            json_lst.append(json.loads(handler.text))
        else:
            json_lst.append('')

    for idx, row in tqdm(dedup_table.iterrows()):
        for i, j in zip(json_lst, urls):

            try:
                evi_handler = str(i['results'][0]['proteinExistence']).split(':')
            except (IndexError, KeyError):
                evi_handler = ['', '']

            dedup_table.loc[idx, 'Protein_existence_evidence'] = evi_handler[1]
            dedup_table.loc[idx, 'Protein_existence_score'] = evi_handler[0]
            dedup_table.loc[idx, 'url_link'] = j

    dedup_table = dedup_table.iloc[:, [0, 1, 12, 11, 2, 3,
                                       4, 5, 6, 7, 8, 9, 10, 13]]

    dedup_table.to_csv('output_with_evidence.tsv', index=False,
                       sep='\t')


if __name__ == '__main__':
    logo()
    if args.check_modules is not None:
        module_checker()
    else:
        if args.make_query is not None:
            if args.make_query == 5:
                query(table=args.make_query[0],
                      name_idx=args.make_query[1],
                      seq_idx=args.make_query[2],
                      output_name=args.make_query[3],
                      delim=args.make_query[5])
            else:
                query(table=args.make_query[0],
                      name_idx=args.make_query[1],
                      seq_idx=args.make_query[2],
                      output_name=args.make_query[3])

        if args.download_TCDB is not None:
            download_db()

        if args.make_blast_db is not None:
            make_db()

            subprocess.run(['./blastp_db_maker.sh'], shell=True)

        if args.blastp is not None:
            make_blastp(query=args.blastp[0],
                        blast_results_name=args.blastp[1])

        if args.make_table is not None:
            remake_table(args.make_table[0],
                         args.make_table[1])

        if args.evidences is not None:
            evidanceer(args.evidences)
