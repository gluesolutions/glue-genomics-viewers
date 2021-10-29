from dataclasses import dataclass
from subprocess import check_call, PIPE, Popen
import os
import logging

import pandas as pd
import pypairix
import tabix
from glue.core import Data

from .subsets import GenomicRangeSubsetState, GenomicMulitRangeSubsetState

from ncls import NCLS #For fast overlaps
from qtpy.QtWidgets import QMessageBox #Ideally QT dependence would go elsewhere


logger = logging.getLogger(__name__)


@dataclass
class GenomeRange:
    chrom: str
    start: int
    end: int

    def __str__(self):
        if self.start <= self.end:
            return f"{self.chrom}:{self.start}-{self.end}"
        else:
            return f"{self.chrom}:{self.end}-{self.start}"


@dataclass
class BedGraph:
    """Utility to aggregate, index and query BedGraph files.

    This will progressively downsample an original dataset, creating a collection
    of aggregation files located in a `.glue_index` folder parallel to the original file.

    The current implementation is restricted to simple, single-attribute files.

    Parameters
    ----------
    path: Path to a local BedGraph file
    downsample_factor: How much each successive aggregation downsamples the previous pass.
    depth: How many passes of downsampling to create and index.
    """
    path: str
    downsample_factor = 10
    depth = 5

    def index(self):
        """
        Aggregate and downsample raw data.

        This operation is currently very slow (10+ minutes / GB)
        """
        os.makedirs(os.path.join(os.path.dirname(self.path), '.glue_index'), exist_ok=True)
        outpaths = [self._level_path(i) for i in range(self.depth)]
        if all(os.path.exists(p + '.bgz.tbi') for p in outpaths):
            logger.debug("Already indexed")
            return

        do_index = show_index_dialog(filename = os.path.basename(self.path))
        print(do_index)
        if do_index == QMessageBox.Ok:
            pass
        else:
            raise FileNotFoundError


        df = pd.read_csv(self.path, delimiter='\t', names=['chr', 'stop', 'start', 'value'])
        check_call(f"bgzip -c {self.path} > {self.path}.bgz", shell=True)
        check_call(['tabix', '-p', 'bed', self.path + '.bgz'])

        downsample_factors = [self.downsample_factor ** (i + 1) for i in range(self.depth)]
        for path, factor in zip(outpaths, downsample_factors):
            decimated = pd.DataFrame(self.decimate((rec for _, rec in df.iterrows()), factor),
                                     columns=['chrom', 'start', 'stop', 'value'])
            decimated.to_csv(path, sep='\t', index=False, header=False)
            df = decimated

            check_call(f"bgzip -c {path} > {path}.bgz", shell=True)
            check_call(['tabix', '-p', 'bed', path + ".bgz"])

    @staticmethod
    def _tabix_query(path, gr: GenomeRange):
        query = str(gr)
        logger.debug("%s %s", path, query)
        try:
            records = tabix.open(path).querys(query)
        except tabix.TabixError:
            return []

        return list(records)

    def query(self, gr: GenomeRange, samples: int = 1000):

        resolution = (gr.end - gr.start) / (samples * 2)
        level = max((i for i in range(self.depth) if self.downsample_factor ** (i + 1) <= resolution),
                    default=None)
        path = self._level_path(level) + '.bgz'
        return pd.DataFrame(self._tabix_query(path, gr), columns=['chrom', 'start', 'stop', 'value']).apply(
            pd.to_numeric, errors='ignore')

    def _level_path(self, level):
        a, b = os.path.split(self.path)

        if level is None:
            return os.path.join(a, '.glue_index', b)

        return os.path.join(a, '.glue_index', '%s.dec_%i_%i' % (b, self.downsample_factor, level))

    @staticmethod
    def decimate(intervals, step):
        try:
            chrom, start, stop, value = next(intervals)
        except StopIteration:
            return

        # TODO: detect non-sorted input

        # TODO: chromosome jumping
        lo, hi = 0, step
        current_sum = 0
        current_max = float('-inf')

        for chrom, start, stop, value in intervals:

            # Start of a new interval
            if lo is None or start >= hi:
                if current_sum > 0:
                    # yield chrom, lo, hi, current_sum / (hi - lo)
                    yield chrom, lo, hi, current_max

                if (stop - start) > step:
                    # spans multiple steps
                    yield chrom, start, stop, value
                    lo = None
                    current_sum = 0
                    current_max = float('-inf')
                else:
                    lo = start
                    hi = lo + step
                    current_sum = value * (stop - start)
                    current_max = max(current_max, value)

            # Straddles an interval
            elif stop > hi:
                current_sum += value * (stop - start)
                # yield chrom, lo, stop, current_sum / (stop - lo)
                yield chrom, lo, stop, max(current_max, value)
                lo = None
                current_sum = 0
                current_max = float('-inf')

            # Contained in an interval
            else:
                current_sum += value * (stop - start)
                current_max = max(current_max, value)

        if current_sum > 0:
            # yield chrom, lo, hi, current_sum / (hi - lo)
            yield chrom, lo, hi, current_max / (hi - lo)


def show_index_dialog(filename):

    msg = QMessageBox()
    msg.setIcon(QMessageBox.Information)
 
    msg.setText(f"The dataset {filename} needs to be indexed by glue before it can be visualized in the Genome Track Viewer.")
    msg.setInformativeText("This process is quite slow (~10min/Gb) but only needs to be done once (indexed files are written to disk in a hidden .glue_index/ directory). Do you want to continue? Glue will be unresponsive during the indexing.")
    msg.setWindowTitle("Do you want to index this file?")
    msg.setStandardButtons(QMessageBox.Ok | QMessageBox.Cancel)
    retval = msg.exec_()
    return retval

@dataclass
class BedPe:
    """
    Utility to aggregate, downsample, and query BedPe files.

    Currently restricted to single-attribute datasets.
    """
    path: str
    downsample_factor = 10
    depth = 7

    def index(self):
        
        os.makedirs(os.path.join(os.path.dirname(self.path), '.glue_index'), exist_ok=True)
        outpaths = [self._level_path(i) for i in range(self.depth)]
        if all(os.path.exists(p + '.bgz.px2') for p in outpaths):
            logger.debug("Already indexed")
            return

        do_index = show_index_dialog(filename = os.path.basename(self.path))
        print(do_index)
        if do_index == QMessageBox.Ok:
            pass
        else:
            raise FileNotFoundError
        
        df = pd.read_csv(self.path, delimiter='\t',
                         names=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'value'])
        df = df.sort_values(['chrom1', 'chrom2', 'start1', 'start2'])

        assert (df.start2 >= df.end1).all()  # all loops are sorted
        assert (df.chrom1 == df.chrom2).all()  # all loops are intra-chromosomal

        check_call(f"sort -k1,1 -k4,4 -k2,2n -k5,5n {self.path} | bgzip  > {self.path}.bgz", shell=True)
        check_call(
            ['pairix', '-f', '-s', '1', '-d', '4', '-b', '2', '-e', '3', '-u', '5', '-v', '6', self.path + '.bgz'])

        downsample_factors = [self.downsample_factor ** (i + 1) for i in range(self.depth)]
        for path, factor in zip(outpaths, downsample_factors):
            df = self.decimate_loops(df, factor)
            df.to_csv(path, sep='\t', index=False, header=False)

            check_call(f"bgzip -c {path} > {path}.bgz", shell=True)
            logger.info("pairix", path)
            check_call(
                ['pairix', '-f', '-s', '1', '-d', '4', '-b', '2', '-e', '3', '-u', '5', '-v', '6', path + '.bgz'])

    @staticmethod
    def _pairix_query(path, gr: GenomeRange, required_endpoints, verbose):
        query = str(gr)

        if required_endpoints == 'either':
            query = [f'{query}|*', f'*|{query}']
        else:
            query = [query]

        if verbose:
            logger.info("%s %s", path, query)

        return list(pypairix.open(path).querys(*query))

    def query(self, gr: GenomeRange, target=100, required_endpoints='both', level=None, verbose=True):
        last = None

        if level is not None:
            path = self._level_path(level) + '.bgz'
            return pd.DataFrame(self._pairix_query(path, gr, required_endpoints, verbose),
                                columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'value']).apply(
                pd.to_numeric, errors='ignore')

        for level in range(self.depth)[::-1]:
            path = self._level_path(level) + '.bgz'
            df = pd.DataFrame(self._pairix_query(path, gr, required_endpoints, verbose),
                              columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'value']).apply(
                pd.to_numeric, errors='ignore')
            if len(df) > target:
                return last if last is not None else df
            last = df

        return df

    def _level_path(self, level):

        a, b = os.path.split(self.path)
        if level is None:
            return os.path.join(a, '.glue_index', b)

        return os.path.join(a, '.glue_index', '%s.dec_%i_%i' % (b, self.downsample_factor, level))

    @staticmethod
    def decimate_loops(df, resolution):
        extent = (df.end2 - df.start1)

        # Remove loops with extents below the resolution
        df = df[extent >= resolution]

        # Merge loops that overlap given the resolution
        a = ((df.end1 + df.start1) // 2) // resolution
        b = ((df.end2 + df.start2) // 2) // resolution
        df = df.groupby([a, b]).apply(lambda x: x.head(1).assign(value=x.value.sum()))

        return df.sort_values(['chrom1', 'chrom2', 'start1', 'start2'])


class GenomicData(Data):
    """Base Data class for wrap Genomic files in Glue."""
    engine_cls = None

    def __init__(self, path, **kwargs):

        kwargs.setdefault('label', os.path.splitext(os.path.split(path)[-1])[0])
        # The task of making a true Glue data object entirely divorced from numpy arrays
        # is involved. As a workaround, we give this dataset a single dummy attribute
        # so that it has a placeholder shape, dimensionality, etc. We largely ignore
        # this information in purpose-built data viewers, but these datasets will have
        # strange behavior if attempting to use standard glue functionality.
        # TODO: Fix this hack
        kwargs['_'] = [0]
        super().__init__(**kwargs)

        self.path = path
        self.engine = self.engine_cls(self.path)

    def profile(self, chr, start, end, subset_state=None, **kwargs):
        raise NotImplementedError


class BedgraphData(GenomicData):
    """Glue Data wrapper for BedGraph files."""

    engine_cls = BedGraph

    def profile(self, chr, start, end, subset_state=None, **kwargs):
        
        query_chrom = f'chr{chr}'
        query_start = int(start)
        query_end   = int(end)
        
        result = self.engine.query(GenomeRange(query_chrom, query_start, query_end))
        if subset_state is None:
            return result

        if isinstance(subset_state, GenomicRangeSubsetState):
            c, s, e = subset_state.chrom, subset_state.start, subset_state.end
            return result.loc[
                result.chrom.eq(c) &
                result.start.ge(s) &
                result.stop.le(e)
            ]
        elif isinstance(subset_state, GenomicMulitRangeSubsetState): 
            mask = result.chrom.eq('junk')
            for c,s,e in zip(subset_state.chroms, subset_state.starts, subset_state.ends):
                if (c != query_chrom) or (s > query_end) or (e < query_start):
                    continue #Avoid making a really long query
                else:
                    mask |= result.chrom.eq(c) & result.start.ge(s) & result.stop.le(e) #This assumes an OR state, but it could be more complicated
            return result.loc[mask]
        else:
            #try:
            #    subset_state.to_index_list()
            #except:
            # TODO: implement more general subset filtering.
            return result.head(0)


class BedPeData(GenomicData):
    """Glue Data wrapper for BedPe files."""
    engine_cls = BedPe

    def profile(self, chr, start, end, subset_state=None, target=100, required_endpoints='both', **kwargs):
        query_chrom = f'chr{chr}'
        query_start = int(start)
        query_end   = int(end)
        result = self.engine.query(
            GenomeRange(query_chrom, query_start, query_end),
            target=target,
            required_endpoints=required_endpoints,
            verbose=False
        )
        if subset_state is None:
            return result

        if isinstance(subset_state, GenomicRangeSubsetState):
            c, s, e = subset_state.chrom, subset_state.start, subset_state.end
            return result.loc[
                (
                    result.chrom1.eq(c) &
                    result.start1.ge(s) &
                    result.end1.le(e)
                ) |
                (
                    result.chrom2.eq(c) &
                    result.start2.ge(s) &
                    result.end2.le(e)
                )
            ]
        elif isinstance(subset_state, GenomicMulitRangeSubsetState):
            # A bit convoluted but reasonably performant
            ncls1 = NCLS(result.start1,result.end1,result.index.values) #Assumes we're in the right chromosome
            ncls2 = NCLS(result.start2,result.end2,result.index.values)
            
            #We could first filter to remove subset_states outside the range entirely
            subset_indices = pd.Series(list(range(len(subset_state.starts)))) #NCLS requires this, though we could do it in the subset object
            subset_starts = pd.Series(subset_state.starts).astype('int')
            subset_ends = pd.Series(subset_state.ends).astype('int')
            
            l_idxs_1, r_idxs_1 = ncls1.all_overlaps_both(subset_starts.values, subset_ends.values, subset_indices.values)
            l_idxs_2, r_idxs_2 = ncls2.all_overlaps_both(subset_starts.values, subset_ends.values, subset_indices.values)
            both_ends = list(set(r_idxs_1) | set(r_idxs_2))
            return result.loc[both_ends]

        else:
            # TODO: implement more general subset filtering.
            return result.head(0)