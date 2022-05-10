from pace_util import resource_list_with_prefix
from pace_util import LayoutManager
from pace_util import Plate96
from pace_util import Plate384
from pace_util import Tip96
from pace_util import LAYFILE
from pace_util import layout_item
from pace_util import ResourceType
from pace_util import ROBOID

lay_mgr = LayoutManager(LAYFILE)

dilution_plate = layout_item(lay_mgr, Plate96, 'dilution_plate')
drug_a_plate = layout_item(lay_mgr, Plate384, 'drug_a_plate')
drug_b_plate = layout_item(lay_mgr, Plate384, 'drug_b_plate')
drug_source_plate = layout_item(lay_mgr, Plate96, 'drug_source_plate')

reader_tray = layout_item(lay_mgr, Plate96, 'reader_tray_'+ROBOID)

waste = layout_item(lay_mgr, Tip96, 'ht_hw_96washdualchamber2_0001')

dilution_tips = layout_item(lay_mgr, Tip96, 'dilution_tips')
pinning_tips = layout_item(lay_mgr, Tip96, 'pinning_tips')
drug_a_tips = layout_item(lay_mgr, Tip96, 'drug_a_tips')
drug_tips_list = [[(drug_a_tips, idx)] + [None]*7 for idx in range(96)]

def drug_tips_iter():
    while True:
        for i in drug_tips_list:
            yield i
drug_tips_iter = drug_tips_iter()


drug_source_tips = layout_item(lay_mgr, Tip96, 'drug_source_tips')