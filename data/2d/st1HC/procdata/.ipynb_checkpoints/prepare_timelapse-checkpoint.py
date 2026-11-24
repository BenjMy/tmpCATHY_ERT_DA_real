from IPython import embed
import duckdb as db
import argparse
#from ertpm.ertio import read_bert
import pandas as pd
import numpy as np
import re
import pygimli as pg

# from https://github.com/Peruz/ERTpm/blob/master/ERTpm/process.py
######################

def read_bert(k_file=None):
    """read bert-type file and return elec and data"""
    with open(k_file) as fid:
        lines = fid.readlines()
    elec_num = int(lines[0])
    data_num = int(lines[elec_num + 2])
    elec_raw = pd.read_csv(k_file, delim_whitespace=True, skiprows=1, nrows=elec_num, header=0)
    elec = elec_raw[elec_raw.columns[:-1]]
    elec.columns = elec_raw.columns[1:]
    data_raw = pd.read_csv(k_file, delim_whitespace=True, skiprows=elec_num + 3, nrows=data_num)
    data = data_raw[data_raw.columns[:-1]]
    data.columns = data_raw.columns[1:]
    return(elec, data)

'''
def to_bert(self, fname, w_rhoa, w_ip, w_err, data_cols, elec_cols, rounding=6):
    if w_rhoa:
        data_cols.extend(['rhoa', 'k'])
    if w_ip:
        data_cols.append('ip')
    if w_err:
        data_cols.append(w_err)
    with open(fname, 'a') as file_handle:
        file_handle.write(str(len(self.elec)) + '\n')
        file_handle.write('# ' + ' '.join(elec_cols) + '\n')
        self.elec[elec_cols].to_csv(file_handle, sep=' ', index=None, header=False)
        # data_wrt = self.data[self.data.valid == 1][data_cols]
        data_wrt = self.data.loc[self.data['valid'] == 1, data_cols]
        file_handle.write(str(len(data_wrt)) + '\n')
        file_handle.write('# ' + ' '.join(data_cols) + '\n')
        data_wrt.to_csv(file_handle, sep=' ', index=None, header=False, float_format='%g')
'''
def to_bert(fname, w_rhoa, w_ip, w_err, data_cols, elec_cols, rounding=6):
    if w_rhoa:
        data_cols.extend(['rhoa', 'k'])
    if w_ip:
        data_cols.append('ip')
    if w_err:
        data_cols.append(w_err)
    with open(fname, 'a') as file_handle:
        file_handle.write(str(len(self.elec)) + '\n')
        file_handle.write('# ' + ' '.join(elec_cols) + '\n')
        self.elec[elec_cols].to_csv(file_handle, sep=' ', index=None, header=False)
        # data_wrt = self.data[self.data.valid == 1][data_cols]
        data_wrt = self.data.loc[self.data['valid'] == 1, data_cols]
        file_handle.write(str(len(data_wrt)) + '\n')
        file_handle.write('# ' + ' '.join(data_cols) + '\n')
        data_wrt.to_csv(file_handle, sep=' ', index=None, header=False, float_format='%g')

######################
        
def general_get_cmd():
    parse = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    main = parse.add_argument_group("main")
    main.add_argument("-fname", type=str, help="data file to process", nargs="+")
    main.add_argument("-dbname", type=str, default='timelapse.db', help="file name for the database")
    main.add_argument("-v", action="store_true", default=False, help="verbose")
    args = parse.parse_args()
    if isinstance(args.fname, str):
        args["fname"] = [args.fname]
    return args


def process(args):

    dbfile = args['dbname']
    conn = db.connect(database=dbfile, read_only=False)

    nns = []

    for i, f in enumerate(args["fname"]):
        nn = f[0: 6]
        nns.append(nn)
        print(f + ' as table ' + nn)

        ds = read_bert(f)
        #newdata = ds.data
        #newelec = ds.elec
        newdata = ds[1]
        newelec = ds[0]

        if args['v']:
            print(newdata)
            print(newelec)

        if i == 0:
            conn.execute(f"""create or replace table Data as (select a, b, m, n, r, k, rhoa, err, '{nn}' as ds from newdata);""")
            conn.commit()
        else:
            conn.execute(f"""insert into Data select a, b, m, n, r, k, rhoa, err, '{nn}' as ds from newdata;""")

    conn.execute("""create or replace view MeanData as (
        select a, b, m, n, count(r) as cnt, mean(r) as r, mean(k) as k, mean(rhoa) as rhoa, mean(err) as err from Data group by a, b, m, n
    );""")

    conn.execute("""create or replace view FullMeanData as (
        select names.ds, a, b, m, n, r, k, rhoa, err, cnt from MeanData cross join (select distinct ds from Data) as names order by ds, a, b, m, n
    );""")

    conn.execute("""create or replace view join_Data_FullMeanData as (
        select
        fmd.ds as ds,
        fmd.a as a,
        fmd.b as b,
        fmd.m as m,
        fmd.n as n,
        Data.ds,
        Data.a,
        Data.b,
        Data.m,
        Data.n,
        Data.r,
        Data.k,
        Data.rhoa as rhoa,
        Data.err as err,
        fmd.r as fmd_r,
        fmd.k as fmd_k,
        fmd.rhoa as fmd_rhoa,
        fmd.err as fmd_err,
        fmd.cnt as fmd_cnt,
        from Data full join FullMeanData as fmd
        using (ds, a, b, m, n)
        order by ds, a, b, m, n
    );""")

    conn.execute("""create or replace table FilledData as (
        select
        ds, a, b, m, n,
        case when r is null then fmd_r else r end as r,
        case when k is null then fmd_k else k end as k,
        case when rhoa is null then fmd_rhoa else rhoa end as rhoa,
        case when err is null then 999.99 else err end as err
        from join_Data_FullMeanData
        where fmd_cnt >= 3
        order by ds, a, b, m, n
    );""")

    #embed()

    for nn, f in zip(nns, args['fname']):
        print(nn, f)
        ds = read_bert(f)
        data_timelapse = conn.execute(f"""
            select a, b, m, n, r, k, rhoa, err
            from filledData
            where ds = '{nn}'
            order by a, b, m, n;
        """).df()
        # print(data_timelapse)
        #data_timelapse = (ds[1]) # non serve
        #ds.data = data_timelapse
        #ds.to_bert(f.replace('.', '.tl'), w_err=True, w_ip=False, w_rhoa=True, w_k=True)
        scheme=pg.DataContainerERT()
        scheme.setSensorPositions(newelec.to_numpy())
        abmn = data_timelapse[["a", "b", "m", "n"]].to_numpy()
        for i, elec in enumerate("abmn"):
            scheme[elec] = abmn[:, i] - 1
        #print(data_timelapse.rhoa)
        #print(data_timelapse.a)
        scheme['rhoa']=data_timelapse["rhoa"]
        scheme['r']=data_timelapse["r"]
        scheme['k']=data_timelapse["k"]
        scheme['err']=data_timelapse["err"]
        scheme.save(f.replace('dat','tldat'))
        print("fullTL file prepare: DONE")
        print(scheme)

if __name__ == '__main__':
    args = vars(general_get_cmd())
    process(args)
