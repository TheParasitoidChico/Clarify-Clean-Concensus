#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec  1 14:07:52 2022
BOLDscraperII library
do not run
contains use case
@author: christian
"""
###############################################################################
runfile('/home/christian/Desktop/COMP/TACHINIDAE_ORD/Host_PD_covariates/TachMatrix.py', wdir='/home/christian/Desktop/COMP/TACHINIDAE_ORD/Host_PD_covariates')
runfile('/home/christian/Desktop/COMP/TachOrd.py', wdir='/home/christian/Desktop/COMP')
###############################################################################
import numpy as np
from Bio import pairwise2
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio import Entrez
import subprocess as sbp
import sqlite3
import os
###############################################################################
############use_case###########################################################
###############################################################################
data = ['Compsilura concinnata', 'Archytas apicifer']
###############################################################################
def ocsv(_path, delim = ',', lower = True):
    env = list()
    if lower == True:
        with open(_path, 'r') as _p:
            for i in _p:
                i = i.strip('\n')
                i = i.split(delim)
                i = [e.lower() for e in i]
                i = [e.strip() for e in i]
                #i = [e.strip('\'') for e in i]
                env.append(i) 
    else:
        with open(_path, 'r') as _p:
            for i in _p:
                i = i.strip('\n')
                i = i.split(delim)
                i = [e.strip() for e in i]
                #i = [e.strip('\'') for e in i]
                env.append(i) 
    return env
###############################################################################
def wcsv(fname, _data, write_type = 'w'):
    with open(fname, write_type) as _f_:
        for _l in _data:
            for l in _l:
                _f_.write(f'{l},')
            _f_.write(f'\n')
###############################################################################
def scrubBOLD(searchterm):
    import requests
    front = 'http://v3.boldsystems.org/index.php/API_Public/sequence?taxon='
    send = front + searchterm
    call = requests.get(send)
    call = call.content
    call = str(call)
    call = call.replace('\\r\\n', '\n')
    call = call.split('\n')
    call = dict(zip(call[::2], call[1::2]))
    return call
###############################################################################
def query_sequence(name, SQL = True, check = 'Arthropoda'):
    bold_query = scrubBOLD(name)
    omit = list()
    for _k, _v in bold_query.items():
        _ksplit = _k.split('|')
        if len(_ksplit) == 4:
            _meta = entrez_taxonomy(entrez_seq_data(_ksplit[3]))
            if check in _meta:
                if SQL == True:
                    SQLsequence([_ksplit[1],_k,_v])
                else:
                    pass
            elif check not in _meta:
                omit.append(_k)#           if !=check, removed.
        elif len(_ksplit) != 4:
            print(f'no accession ID - for now, removed:{_k}')
            omit.append(_k)
    for omitted_k in omit:
        bold_query.pop(omitted_k)
    return bold_query
###############################################################################
def entrez_seq_data(genbank_accession_id, db = 'nucleotide', email = 'christianconnors@nevada.unr.edu'): # please, use your own email
    Entrez.email = email #required for API compliance
    handle = Entrez.efetch(db=db, id = genbank_accession_id, rettype = 'gb', retmode = 'text')
    ohandle = SeqIO.read(handle, 'genbank')
    return ohandle
###############################################################################
def entrez_taxonomy(seqio_obj, annotation_field = 'taxonomy'):
    extract = seqio_obj.annotations.get(annotation_field)
    return extract
###############################################################################
def SQLsequence(data, db = 'calls.db'):
    connection = sqlite3.connect(db)
    cursor = connection.cursor()
    SQL_kv_update = "INSERT INTO call (taxon,gene,id,accessionid,sequence) VALUES(?,?,?,?,?)"
    TO_SQL = list()
    if isinstance(data[0], list):
        for each in data:
            _tuple = (*each,)
            TO_SQL.append(_tuple)
        cursor.executemany(SQL_kv_update, TO_SQL)
        connection.commit()
    elif isinstance(data[0], str):
        _tuple = (*data,)
        TO_SQL.append(_tuple)
        cursor.executemany(SQL_kv_update, TO_SQL)
        connection.commit()
    elif isinstance(data, tuple):
        TO_SQL.append(data)
        cursor.executemany(SQL_kv_update, TO_SQL)
        connection.commit()
    else:
        print(f'Error: input != tuple, list, or simple array')
###############################################################################
def wfa(fa, outpath):
    for key, val in fa.items():
        key = key.replace(':', '')
        key = key.replace('.', '')
        key = key.replace(',', '')
        key = key.replace('b\'>', '>')
        key = dewhite(key)
        print(key)
        with open(outpath, 'a') as OUT:
            OUT.write(f'{key}\n{val}\n')
###############################################################################
def mafft(path_fa, out, in_house = False):
    print('mafft\nmultiple sequence alignment with fast fourier transformation')
    mafft_command = '"/usr/bin/mafft"  --auto --clustalout --reorder {fa} > {out}'.format(fa=path_fa, out=out)
    sbp.call(mafft_command, shell=True)
###############################################################################
def oaln(path_aln, aln_format = 'clustal'):
    aln = AlignIO.read(path_aln, aln_format)
    return aln
###############################################################################
# need to ignore (-) unless there are only missing characters.
# it ... should do that. it doesn't seem to for some reason
def conseq(aln, test = list()):
    untransposed = list()
    concensus = list()
    for sequence in aln:
        untransposed.append(sequence.seq)
    transposed = np.array(untransposed).T.tolist()
    for i in transposed:
        max_nucleotide = max(set(i), key = i.count)
        if max_nucleotide != '-':
            concensus.append(max_nucleotide)
        elif max_nucleotide == '-':
            i_omitted_missing = [n for n in i if n != '-']
            if len(i_omitted_missing) == 0:
                concensus.append('-')
            elif len(i_omitted_missing) >= 1:
                max_nucleo_omit = max(set(i_omitted_missing), key = i_omitted_missing.count)
                concensus.append(max_nucleo_omit)
    con_seq = str()
    con_seq = con_seq.join(concensus)
    return con_seq
###############################################################################
def query_sequence(name, SQL = True, check = 'Arthropoda'):
    print('QUERY SEQUENCE')
    bold_query = scrubBOLD(name)
    omit = list()
    for _k, _v in bold_query.items():
        _ksplit = _k.split('|')
        if len(_ksplit) == 4:
            _meta = entrez_taxonomy(entrez_seq_data(_ksplit[3]))
            if check in _meta:
                if SQL == True:
                    SQLsequence([_ksplit[1],_k,_v])
                else:
                    pass
            elif check not in _meta:
                omit.append(_k)#           if !=check, removed.
        elif len(_ksplit) != 4:
            print(f'no accession ID - for now, removed:{_k}')
            omit.append(_k)
    for omitted_k in omit:
        bold_query.pop(omitted_k)
    return bold_query
###############################################################################
def dict_concensus(dict_fa, keep = False):
    print('DICT CONCENSUS')
    removed_files = ['temp.fa', 'tempout']
    if keep == False:
        if 'temp' not in os.listdir():
            wfa(dict_fa, 'temp.fa')
            mafft('temp.fa', 'tempout')
            aligned = oaln('tempout')
            alig_seq = conseq(aligned)
            for removed in removed_files:
                os.remove(removed)
            return alig_seq
        else:
            print('Error: temp already exists/n...removing temp files...')
            for removed in removed_files:
                os.remove(removed)
###############################################################################
def if_not_in(item, itemlist):
    if item not in itemlist:
        itemlist.append(item)
###############################################################################
def if_already_db(obj, db, element = 0, return_row = False):
    print('IF ALREADY IN DB')
    in_db = list()
    element_db = list()
    connection = sqlite3.connect(db)
    _cursor = connection.execute("SELECT * from call WHERE taxon = ?", tuple([obj]))
    [in_db.append(i) for i in _cursor]
    if return_row == True:
        return in_db
    elif return_row == False:
        for _row_returned in in_db:
            if_not_in(_row_returned[element], element_db)
        if obj not in element_db:
            return False
        elif obj in element_db:
            return True
###############################################################################
def CCC(data, db): # output should be dictionary to deal with multiple genes per taxa
    if isinstance(data, str):
        cumulative = list()
        if if_already_db(data, db):
            taxon = data
            connection = sqlite3.connect(db)
            _cursor = connection.cursor()
            db_genes = [c[0] for c in _cursor.execute("SELECT gene FROM call")]
            align_genes = list()
            [if_not_in(gene, align_genes) for gene in db_genes]
            for gene in align_genes:
                c = call_and_concensus(taxon, gene, db)
                cumulative.append(c)
            return cumulative
        else:
            submissions = create_SQL_row(query_sequence(data, False, 'Arthropoda'))
            for taxon in submissions[0]:
                for gene in submissions[1]:
                    c = call_and_concensus(taxon, gene, db)
                    cumulative.append(c)
            return cumulative
###############################################################################
def create_SQL_row(BOLDdict):
    print('CREATE SQL ROW')
    taxon_collector = list()
    gene_collector = list()
    for bold_key, seq_val in BOLDdict.items():
        bold_key_split = bold_key.split('|')
        taxon = bold_key_split[1]
        if_not_in(taxon, taxon_collector)
        gene = bold_key_split[2]
        if_not_in(gene, gene_collector)
        accessionid = bold_key_split[3]
        _out_tuple = tuple([taxon,gene,bold_key,accessionid,seq_val])
        print(_out_tuple)
        SQLsequence(_out_tuple)
    return [taxon_collector, gene_collector]
###############################################################################
def call_and_concensus(taxon, gene, db):
    print('CALL AND CONCENSUS')
    taxon_gene_dict = dict()
    connection = sqlite3.connect(db)
    _cursor = connection.execute("SELECT * FROM call WHERE taxon = ? AND gene = ?", tuple([taxon,gene]))
    for I in _cursor:
        sub_dict = {I[2]:I[4]}
        taxon_gene_dict.update(sub_dict)

    if len(taxon_gene_dict) == 1:
        print(taxon_gene_dict)
        _out = [taxon, gene]
        for key, val in taxon_gene_dict.items():
            seq = taxon_gene_dict.get(key)
            seq = seq.lower()
            _out.append(seq)
        return tuple(_out)

    elif len(taxon_gene_dict) > 1:
        concensus_seq = dict_concensus(taxon_gene_dict)
        taxon_concensus = tuple([taxon,gene,concensus_seq])
        return taxon_concensus
    
    elif len(taxon_gene_dict) == 0:
        pass
###############################################################################
def commit_conseq(_conseqtuple, db):
    #SQL_kv_update = "INSERT INTO consensus_sequences (taxon,gene,concensus) VALUES(?,?,?)"
    for _conseq in _conseqtuple:
        if isinstance(_conseq, tuple):
            print(_conseq)
            connection = sqlite3.connect(db)
            cursor = connection.cursor()
            cursor.execute("INSERT INTO consensus_sequences (taxon,gene,concensus) VALUES(?,?,?)", _conseq)
            connection.commit()
###############################################################################
def gbify(searchterm, *args):
    