import sys, os, re, copy
import numpy as np, pandas as pd, tables as tb
tb.parameters.MAX_GROUP_WIDTH = 51200
# tb.parameters.NODE_CACHE_SLOTS = -51200
# tb.parameters.METADATA_CACHE_SIZE = 1048576 * 100000
# tb.parameters.CHUNK_CACHE_SIZE = 2097152 * 100000
# tb.parameters.CHUNK_CACHE_NELMTS = 521
import scipy.stats as st

class Environment:
    def __init__(self):
        self.float = np.float32
        self.duplicate_tag = '_duplicated_'
        self.common_suffix = '.h5'

env = Environment()

class TBData(dict):
    def __init__(self, data, name, msg = None, root = '/', complib = 'bzip2'):
        '''bzip2 may not be compatible with other hdf5 applications; but zlib is fine'''
        self.__root = root.strip('/')
        self.__group = name
        self.__msg = msg
        try:
            if type(data) is dict:
                self.update(data)
            elif type(data) is str:
                # is file name
                self.__load(tb.open_file(data))
            else:
                # is file stream
                self.__load(data)
        except tb.exceptions.NoSuchNodeError:
            raise ValueError('Cannot find dataset {}!'.format(name))
        self.tb_filters = tb.Filters(complevel = 9, complib=complib)

    def sink(self, filename):
        with tb.open_file(filename, 'a') as f:
            if self.__root:
                try:
                    f.create_group("/", self.__root)
                except:
                    pass
            try:
                # there is existing data -- have to merge with current data
                # have to do this because the input file lines are not grouped by gene names!!
                # use try ... except to hopefully faster than if ... else
                # e.g., if not f.__contains__('/{}'.format(self.__group)) ... else ...
                for element in f.list_nodes('/{}/{}'.format(self.__root, self.__group)):
                    if element.name != 'colnames':
                        self[element.name] = np.concatenate((element[:], self[element.name]))
            except tb.exceptions.NoSuchNodeError:
                f.create_group("/" + self.__root, self.__group,
                               self.__msg if self.__msg else self.__group)
            for key in self:
                self.__store_array(key, f)
            f.flush()

    def dump(self, table, output = False):
        if output:
            pd.DataFrame({self['colnames'][i] : self[table][:,i] for i in range(len(self['colnames']))}, index = self['rownames']).to_csv(sys.stdout, na_rep = 'NA')
            return None
        else:
            return pd.DataFrame({self['colnames'][i] : self[table][:,i] for i in range(len(self['colnames']))}, index = self['rownames'])

    def __load(self, fstream):
        try:
            for element in fstream.list_nodes('/{}/{}'.format(self.__root, self.__group)):
                self[element.name] = element[:]
            fstream.close()
        except:
            fstream.close()
            raise

    def __roll_back(self, group, name):
        try:
            n = getattr(group, name)
            n._f_remove()
        except AttributeError:
            pass

    def __store_array(self, name, fstream):
        if self.__root:
            element = getattr(getattr(fstream.root, self.__root), self.__group)
        else:
            element = getattr(fstream.root, self.__group)
        arr = self[name]
        if type(arr) is list:
            arr = np.array(arr)
        self.__roll_back(element, name)
        #
        if arr.shape != (0,):
            ds = fstream.create_carray(element, name, tb.Atom.from_dtype(arr.dtype), arr.shape,
                                       filters = self.tb_filters)
            ds[:] = arr

def get_tb_grps(filenames, group_name = None):
    if len(filenames) == 0:
        return []
    if isinstance(filenames, str):
        filenames = [filenames]
    names = set()
    for filename in filenames:
        with tb.open_file(filename) as f:
            names.update([node._v_name for node in (f.root if group_name is None else getattr(f.root, '{}'.format(group_name)))])
    return sorted(names)

def get_se(bhat, pval):
    if bhat < 0:
        z = st.norm.ppf(pval / 2)
    else:
        z = st.norm.ppf(1 - pval / 2)
    if z != 0:
        se = bhat/z
    else:
        se = np.nan
    return se

class SSData:
    def __init__(self, header = False):
        self.data = {'buffer':{'data':[], 'rownames':[]}, 'output':{}}
        self.header = header
        self.previous_name = self.current_name = None
        self.min_line_len = max([3, 4, 5])

    def parse(self, line, ensg_version = 0, include = set()):
        # input line is snp, gene, beta, t, pval and optionally se(bhat)
        if not line:
            self.__reset()
            self.current_name = None
            return 1
        if isinstance(line, bytes):
            line = line.decode()
        line = line.strip().split()
        if self.header:
            # header line has to be accounted for
            self.header = False
            return 0
        #
        line[0] = line[0].strip('"')
        if not line[0] in include and len(include)>0:
            return -1
        # the line is not long enough ... there is format issues
        if len(line) < self.min_line_len:
            return -1
        #
        if ensg_version == 0:
            line[0] = line[0].rsplit('.',1)[0]
        if self.previous_name is None:
            self.previous_name = line[0]
        self.current_name = line[0]
        if self.current_name != self.previous_name:
            self.__reset()
        if 5 == 0:
            # the two numbers are bhat and p-value, we need to compute se(bhat)
            self.data['buffer']['data'].append([line[2], str(get_se(float(line[2]), float(line[3]))), line[3]])
        else:
            self.data['buffer']['data'].append([line[2], line[3], line[4]])
        self.data['buffer']['rownames'].append(self.__format_variant_id(line[1]))
        return 0

    @staticmethod
    def __format_variant_id(value):
        value = value.strip('"').lstrip('chr').split('_')
        if len(value) > 4:
            # keep it chr, pos, ref, alt
            value = value[:4]
        return '_'.join(value)

    def __reset(self):
        self.data['buffer']['data'] = np.array(self.data['buffer']['data'], dtype = env.float)
        self.data['buffer']['rownames'] = np.array(self.data['buffer']['rownames'])
        self.data['buffer']['colnames'] = np.array(['beta','se','pval'])
        self.data['output'] = copy.deepcopy(self.data['buffer'])
        self.data['buffer'] = {'data':[], 'rownames':[]}

    def dump(self):
        return self.data['output']

class DataMerger(TBData):
    def __init__(self, files, name, msg = None):
        TBData.__init__(self, {}, name, msg, complib = "zlib")
        self.files = sorted(files)
        self.__group = name

    def merge(self):
        data = {}
        one_snp = None
        failure_ct = 0
        # Collect data
        for item in self.files:
            tissue = re.sub(r'{}$'.format(env.common_suffix), '', os.path.basename(item))
            try:
                data[tissue] = TBData(item, self.__group)
                if one_snp is None: one_snp = data[tissue]['rownames'][0]
            except ValueError:
                data[tissue] = {'data' : np.array([[np.nan, np.nan, np.nan]]), 'rownames': None}
                failure_ct += 1
            # Fix row name
            # Because in GTEx data file there are duplicated gene-snp pairs having different sumstats!!
            if data[tissue]['rownames'] is not None:
                data[tissue]['rownames'] = self.__dedup(data[tissue]['rownames'], item)
        if failure_ct == len(self.files):
            return 1
        # Merge data
        for idx, item in enumerate(['beta','se','pval']):
            self[item] = pd.concat([pd.DataFrame(
                {tissue : data[tissue]['data'][:,idx]},
                index = data[tissue]['rownames'] if data[tissue]['rownames'] is not None else [one_snp]
                ) for tissue in sorted(data.keys())], sort=True, axis = 1)
            if 'rownames' not in self:
                self['rownames'] = np.array(self[item].index, dtype = str)
            if 'colnames' not in self:
                self['colnames'] = np.array(self[item].columns.values.tolist(), dtype = str)
            self[item] = np.array(self[item].values, dtype = env.float)
        # np.savetxt(sys.stdout, self['pval'], fmt='%10.5f')
        # print(self['rownames'])
        # print(self['colnames'])
        return 0

    def __dedup(self, seq, filename):
        seen = {}
        dups = set()
        def __is_seen(x, seen):
            if x not in seen:
                seen[x] = 0
                return 0
            else:
                seen[x] += 1
                dups.add(x)
                return 1
        # Tag them
        obs = [x if not __is_seen(x, seen) else '%s%s%s' % (x, env.duplicate_tag, seen[x]) for x in seq]
        # Log them
        if len(dups):
            filename = os.path.splitext(filename)[0]
            with open(filename + '.error', 'a') as f:
                for item in dups:
                    f.write('{}:{} appeared {} times in {}\n'.\
                            format(self.__group, item, seen[item] + 1, filename))
        return obs

def get_gs_pairs(data, name, num = (1, 9), method = 'equal_space'):
    '''choose gene-snp pairs from data, controlled by num = (a, b)
    for the best gene-snp pair (a = 0 or 1), and b other random 
    gene-snp pairs'''

    def random_sample(x, k):
        if k < 0:
            return sorted(x)
        return sorted(set(np.random.choice(x, min(len(x), k))))

    def equal_sample(x, k):
        if len(x) < k or k < 0:
            return sorted(x)
        f = lambda m, n: [i*n//m + n//(2*m) for i in range(m)]
        return sorted(set([x[i] for i in f(k, len(x))]))

    output = {'colnames' : data['colnames']}
    lp = data.dump('pval')
    shat = data.dump('se')
    # pval cannot be nan or zero (in fastqtl output, Shat nan has pval zero)
    lp = lp[np.all(np.isfinite(lp), axis=1) & np.all(np.isfinite(shat), axis=1)]
    #
    if lp.empty and num[0] > 0:
        return None
    # Find strong SNP-gene pair
    lp = -np.log10(lp)
    rowidx = np.where(data['rownames'] == lp.max(axis=1).idxmax())[0][0]
    if num[0] != 0:
        output['strong_rownames'] = ['%s_%s' % (name, data['rownames'][rowidx].decode())]
        output['strong'] = {}
        for k in ['beta', 'pval', 'se']:
            output['strong'][k] = data[k][rowidx, :]
    if num[1] != 0:
        all_randomidxes = [y for y, x in enumerate(data['rownames']) if x in lp.index and y != rowidx]
        sample_randomidxes = random_sample(all_randomidxes, num[1]) if method == 'random' else equal_sample(all_randomidxes, num[1])
        output['random_rownames'] = ['%s_%s' % (name, data['rownames'][x].decode()) for x in sample_randomidxes]
        output['random'] = {}
        for k in ['beta', 'pval', 'se']:
            output['random'][k] = data[k][sample_randomidxes, :]
    if not 'strong' in output and not 'random' in output:
        output = None
    return output

def merge_tmp_h5(output, verbose = 0):
    from glob import glob
    tmpfiles = list(glob(output + "_*.tmp"))
    if os.path.isfile(output+'.h5'):
        os.remove(output+'.h5')
    for item in sorted(tmpfiles):
        for name in get_tb_grps(item):
            cmd = 'h5copy -i {0} -o {2} -s "/{1}" -d "/{1}"'.format(item, name, output+'.h5')
            if verbose:
                print(cmd)
            os.system(cmd)

