"""Parser specification checks; synthetic values are not performance evidence."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
sys.dont_write_bytecode = True
ROOT = Path(__file__).resolve().parents[4]
spec = importlib.util.spec_from_file_location('lane', ROOT/'scripts/lane_verify_parse.py')
lane = importlib.util.module_from_spec(spec)
spec.loader.exec_module(lane)

class ParserTests(unittest.TestCase):
    def record(self, rank):
        return dict(rank=rank, iteration=4, cutoff_ha=2.5, event='attempt', status=1, timed=True,
                    sizes=dict.fromkeys(lane.SIZES, 1), seconds=dict.fromkeys(lane.PARTS, 0.),
                    map_seconds=dict.fromkeys(lane.MAP_PARTS, 0.), h_vectors=7, davidson_iterations=2)
    def parse(self, records):
        with tempfile.TemporaryDirectory(dir=ROOT/'target/scratch') as d:
            f = Path(d)/'results.json'
            f.write_text(json.dumps(dict(phase_profile=dict(ranks=2, verification_lane=dict(schema='verify-profile-1', records=records)))))
            return lane.parse(f)
    def test_critical_rank_uses_closed_rank_total_not_sum_of_maxima(self):
        a,b = self.record(0),self.record(1)
        a['seconds'].update(warm_start=9.,map=1.)
        b['seconds'].update(warm_start=1.,map=8.)
        b['map_seconds']['solve']=7.
        rows, groups = self.parse([a,b])
        self.assertEqual(groups[0]['critical_rank'],0)
        self.assertEqual(groups[0]['max_total_s'],10.)
        self.assertEqual(groups[0]['min_total_s'],9.)
        self.assertEqual(rows[0]['warm_copy_read_write_bytes'],2)
        self.assertIn('Critical rank', lane.render(rows,groups))
    def test_run_total_does_not_sum_per_attempt_maxima(self):
        records = [self.record(0), self.record(1), self.record(0), self.record(1)]
        for r,t in zip(records,[10.,9.,1.,3.]): r['seconds']['map']=t
        for r in records[2:]: r['iteration']=5
        rows,_ = self.parse(records)
        run = lane.aggregate(rows)[0]
        self.assertEqual(run['critical_rank'],1)
        self.assertEqual(run['max_lane_total_s'],12.)

    def test_count_only_and_schedule_are_separate(self):
        a,b = self.record(0),self.record(0)
        a['timed']=False
        b['event']='schedule'
        b['timed']=False
        rows,groups=self.parse([a,b])
        self.assertEqual(len(groups),2)
        self.assertTrue(all(r['total_s']==0 for r in rows))
        b['seconds']['schedule']=0.1
        with self.assertRaises(ValueError): self.parse([b])
    def test_invalid_payloads_reject(self):
        for key,value in [('rank',2),('status',5),('event','other'),('timed',1)]:
            r=self.record(0); r[key]=value
            with self.assertRaises(ValueError): self.parse([r])
        for value in [float('nan'),float('inf'),-1.]:
            r=self.record(0);r['seconds']['map']=value
            with self.assertRaises(ValueError): self.parse([r])
        r=self.record(0);del r['sizes']['grid_x']
        with self.assertRaises(KeyError): self.parse([r])
    def test_empty_records_are_valid(self):
        self.assertEqual(self.parse([]),([],[]))

if __name__=='__main__': unittest.main()
