# Get Start

```
pea = createPea(n, zone='DE ', resolution = 1)
carriers.power

v2g = pea.get(carriers.power, 'Generator')

v2g = pea.getIn('AC')
v2g = pea.getOut('DC')

v2g.df  
v2g.t('p').sum()
v2g.capex()
v2g.opex()
v2g.export()
v2g.import()

```