import pandas as pd
import yaml
pd.options.mode.chained_assignment = None

seasons = {
    'all': [1,2,3,4,5,6,7,8,9,10,11,12],
    'win': [10,11,12,1,2,3],
    'sum': [4,5,6,7,8,9]
           }

#stat austria
inhaut = 8773000000

mundata = {
    'kamptal': {
        'Grafenegg': 31398,
        'Gföhl': 31311,
        'Hadersdorf-Kammern': 31315, #notour
        'Jaidhof': 31319, #notour
        'Langenlois': 31322,
        'Lengenfeld': 31323,
        'St. Leonhard am Hornerwald': 31340, #notour
        'Schönberg am Kamp': 31355,
        'Droß': 31356 #notour
    },
    'nockberge': {
        'Deutsch Griffen': 20503,
        'Bad Kleinkirchheim': 20601,
        'Albeck': 21001,
        'Gnesau': 21004,
        'Reichenau': 21007
        },
    'joglland': {
        'Fischbach': 61708,
        'Miesenbach bei Birkfeld': 61728,
        'Ratten': 61741,
        'Rettenegg': 61743,
        'St. Kathrein am Hauenstein': 61744,
        'Strallegg': 61750,
        'Birkfeld': 61757,
        'St. Jakob im Walde': 62242,
        'Wenigzell': 62262,
        'Vorau': 62278,
        'Waldbach-Mönichwald': 62279
        }
}

#stat austria
inhabitants = {
    'kamptal': [3061,3783,2009,1225,7609,1412,1119,1678,1903,983],
    'nockberge': [908,1723,995,1034,1822],
    'joglland': [1528,696,1122,741,633,1925,4992,1056,1391,4731,1523]
}

# Statistik Austria - Ankünfte und Übernachtungen - Berichtsgemeinden Sommer - Winter 2017/2018
tourism = {
    'stays': {
        'kamptal': {'sum': [3447, 402, 41478, 1610, 8788, 6528],
                    'win': [11423, 787, 20320, 423, 2996, 1283],
                    'all': [3447, 402, 41478, 1610, 8788, 6528, 11423, 787, 20320, 423, 2996, 1283]
                   },
        'nockberge': {'sum': [5003, 319673, 16745,  4702,  91835],
                      'win': [1410, 433398, 18249, 12162, 118617],
                      'all': [5003, 319673, 16745,  4702,  91835, 1410, 433398, 18249, 12162, 118617]
                     },
        'joglland': {
                     'win': [14313,  8020,  4881, 1490, 3604, 1950,  995, 15077, 11259,  3006,  7991],
                     'sum': [20382, 17951, 11616, 2632, 5885, 3890, 5157, 22249, 19053, 10212, 19689],
                     'all': [14313,  8020,  4881, 1490, 3604, 1950,  995, 15077, 11259,  3006,  7991, 20382, 17951, 11616, 2632, 5885, 3890, 5157, 22249, 19053, 10212, 19689]
                    }
    },
    'ratio': {
        'kamptal': {'sum': 1, 'win': 1, 'all': 1},
        'nockberge': {'sum': 1.45, 'win': 2.55, 'all': 2.00},
        'joglland': {'sum': 1, 'win': 1, 'all': 1}
    }
}

#e-control
ldaut = {
    1920:  1.763,
    1925:  2.110,
    1930:  2.367,
    1935:  2.132,
    1940:  3.292,
    1945:  2.785,
    1950:  5.640,
    1955:  9.595,
    1960: 13.315,
    1965: 17.816,
    1970: 23.908,
    1975: 30.275,
    1980: 37.473,
    1985: 41.844,
    1990: 48.529,
    1995: 52.606,
    2000: 58.512,
    2001: 60.347,
    2002: 61.073,
    2003: 63.308,
    2004: 64.894,
    2005: 66.083,
    2006: 67.373,
    2007: 67.883,
    2008: 68.516,
    2009: 65.882,
    2010: 68.931,
    2011: 68.992,
    2012: 69.630,
    2013: 69.934,
    2014: 68.942,
    2015: 69.917,
    2016: 70.740,
    2017: 71.824
          }

# entso-e data in MW on an hourl basis (columns), each row is a day
ldata = pd.read_excel('/home/cmikovits/my-data/Monthly-hourly-load-values_2006-2015.xlsx', skiprows = 3, index_col = 0, sheet_name='MHLV_2006-2015')
atldata = ldata.filter(regex='^AT$', axis=0)
atldata['Load'] = atldata.iloc[:,-24:].sum(axis=1)
monthdata = atldata.Load.groupby([atldata.Month]).sum()
nryears = len(atldata.Year.unique())

l = {}
rel = {}
load = {}

for y in ldaut:
    load[y] = {}

#seasonal TWh/a
for s in seasons:
    l[s] = float('{:.3f}'.format(monthdata.loc[seasons[s]].mean()/nryears*12/1000000))
    rel[s] = float('{:.4f}'.format(l[s]/l['all']))
    for y in ldaut:
        load[y][s] = float('{:.3f}'.format(ldaut[y]*rel[s]))

ly = max(ldaut, key=int)
regions = {}

for r in inhabitants:
    regions[r] = {}
    regions[r]['load'] = {}
    regions[r]['stays'] = {}
    regions[r]['stayrelinh'] = {}
    regions[r]['inh'] = int(sum(inhabitants[r]))

    for s in seasons:
        regions[r]['load'][s] = float('{:.0f}'.format(ldaut[ly]*rel[s]/inhaut*regions[r]['inh']*10**9))
        n = 300
        if s != 'all':
            n = 150
        regions[r]['stays'][s] = int(sum(tourism['stays'][r][s])/n)
        regions[r]['stayrelinh'][s] = float('{:.2f}'.format(regions[r]['stays'][s]/regions[r]['inh']*100))




print('kWh per inhabtiant: ', ldaut[ly]/inhaut*10**12)
print('load: ', yaml.dump(load, allow_unicode=True, default_flow_style=False))
print('ratio win sum: ', rel)
print('regions: ', yaml.dump(regions, allow_unicode=True, default_flow_style=False))

