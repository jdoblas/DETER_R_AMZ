import ee
import math

from lib.sar_ee_utils import toDB,toNatural


def harmonic_detrending(col, band):

    col = col.map(toDB)

    def addVariables(image):
        years = image.date().difference('2016-01-01', 'year')
        return image.addBands(ee.Image(years).rename('t')).addBands(ee.Image.constant(1)).float()

    col_st2 = col.select([band]).map(addVariables)
    print(col_st2.size().getInfo())
    # harmonicIndependents = ee.List(['constant', 't', 'cos', 'sin']) # Won't use trend
    harmonicIndependents = ee.List(['constant', 'cos', 'sin'])
    dependent = band

    def addHarmonizedVars(image):
        timeRadians = image.select('t').multiply(2 * math.pi)
        return image.addBands(timeRadians.cos().rename('cos')).addBands(timeRadians.sin().rename('sin'))

    harmonicsS1 = col_st2.map(addHarmonizedVars)
    harmonicLR = harmonicsS1 \
        .select(harmonicIndependents.add(dependent)) \
        .reduce(ee.Reducer.linearRegression(harmonicIndependents.length(), 1))
    # Turn the array image into a multi-band image of coefficients.
    harmonicLRCoefficients = harmonicLR.select('coefficients') \
        .arrayProject([0]) \
        .arrayFlatten([harmonicIndependents])

    def computeFittedHarmonicNoIntercept(image):
        # fit = image.select(harmonicIndependents)\
        #         .multiply(harmonicLRCoefficients)\
        fit = image.select(['cos', 'sin']) \
            .multiply(harmonicLRCoefficients.select(['cos', 'sin'])) \
            .reduce('sum') \
            .rename('fitted')
        return image.addBands(fit) \
            .addBands(image.select(band).subtract(fit).rename("residual"))

    fittedHarmonic = harmonicsS1.map(computeFittedHarmonicNoIntercept)
    col_stable = fittedHarmonic.select("residual")
    return (col_stable.map(toNatural))
