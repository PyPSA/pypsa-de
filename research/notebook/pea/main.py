import pandas as pd


class Grouper:
    def __init__(self, n, df, config={}):
        self.n = n
        self.df = df
        self.config = config

    def exportP(self, p0="p0", p1="p1"):
        advance = self.t(p0).clip(lower=0)
        revert = self.t(p1).clip(lower=0)
        return pd.concat([advance, revert], axis=1)

    def importP(self, p0="p0", p1="p1"):
        advance = -self.t(p0).clip(upper=0)
        revert = -self.t(p1).clip(upper=0)
        return pd.concat([advance, revert], axis=1)

    def t(self, attr, pAs=None):
        componentMethods = [
            "generators_t",
            "links_t",
            "lines_t",
            "stores_t",
            "storage_units_t",
            "loads_t",
        ]

        resultDf = None
        for method in componentMethods:
            dict = getattr(self.n, method)

            currentAttr = attr
            if (method == "links_t" or method == "lines_t") and (attr == "p") and pAs:
                currentAttr = pAs

            if currentAttr not in dict:
                continue

            df = dict[currentAttr]
            if self.df is None:
                continue
            keys = self.df.index
            valid_keys = [key for key in keys if key in df.columns]
            if not valid_keys:
                continue
            currentDf = df[valid_keys]  # * self.config['resolution']

            if resultDf is None or resultDf.empty:
                resultDf = currentDf
            else:
                # if resultDf is not None:
                #   common_indexes = resultDf.index.intersection(currentDf.index)
                #   if not common_indexes.empty:
                #     raise 'common indexes'

                resultDf = pd.concat([resultDf, currentDf], axis=1)

        return resultDf

    def capex(self):
        return self.df["capital_cost"].sum()

    def p_nom_opt(self):
        return self.df["p_nom_opt"].sum()

    def e_nom_opt(self):
        return self.df["e_nom_opt"].sum()

    def energy(self, p="p"):
        energy = self.t(p).sum(axis=1).sum() * self.config["resolution"]
        return abs(energy)


class Pea:
    def __init__(self, n, config={}):
        config["zone"] = config.get("zone", "DE0 ")
        config["resolution"] = config.get("resolution", 1)

        self.n = n
        self.config = config
        self.zone = config.get("zone", "DE0 ")
        self.resolution = config.get("resolution", 1)

    def get(self, carriers, type="common", components=None, filter=None):
        componentMethods = components
        if components == None and type == "common":
            componentMethods = [
                "generators",
                "links",
                "lines",
                "stores",
                "storage_units",
                "loads",
            ]

        if components == None and type != "common":
            componentMethods = ["links", "lines"]

        # componentMethods = ['stores']
        if isinstance(carriers, str):
            carriers = [carriers]

        resultDf = None

        for method in componentMethods:
            df = getattr(self.n, method)

            zone = df.index.str.startswith(self.zone)

            if type == "import":
                zone = (~df["bus0"].str.startswith(self.zone)) & df[
                    "bus1"
                ].str.startswith(self.zone)

            if type == "export":
                zone = df["bus0"].str.startswith(self.zone) & ~df[
                    "bus1"
                ].str.startswith(self.zone)
            dfIndexs = zone & df["carrier"].isin(carriers)

            if type == "inner":
                zone = df["bus0"].str.startswith(self.zone) & df["bus1"].str.startswith(
                    self.zone
                )
            dfIndexs = zone & df["carrier"].isin(carriers)

            if callable(filter):
                dfIndexs = dfIndexs & filter(df)

            currentDf = df[dfIndexs]

            if len(currentDf) > 0:
                if resultDf is None:
                    resultDf = currentDf
                else:
                    if resultDf is not None:
                        common_indexes = resultDf.index.intersection(currentDf.index)
                        if not common_indexes.empty:
                            raise "common indexes"

                    resultDf = pd.concat([resultDf, currentDf], axis=0)
        group = Grouper(self.n, resultDf, self.config)
        return group
