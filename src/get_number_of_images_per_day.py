import ee, datetime
from get_image_ids import get_image_ids
from src.get_config import get_config
import seaborn as sns
import matplotlib.pyplot as plt

ee.Initialize()
AOI = ee.FeatureCollection('users/detersaree/DETER_SAR_ANCILLARY/lm_bioma_250_IBGE2019_AMZ_simplify500').first().geometry()
config = get_config()

date0 = ee.Date('2020-01-01')
total = []
for i in range(365):
    n = len(get_image_ids(date0.advance(i,'days').format('YYYY-MM-dd').getInfo(), date0.advance(i+1,'days').format('YYYY-MM-dd').getInfo(), config))
    total.append(n)

print (total)
sns.distplot(total, kde=True, rug=False)

plt.show()