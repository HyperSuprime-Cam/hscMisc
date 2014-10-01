from .buildAndCatalog import BuildAndCatalog

FILTERS = "UGRIZ"

class BuildSdss(BuildAndCatalog):
    def __init__(self, *args, **kwargs):
        super(BuildSdss, self).__init__(*args, **kwargs)
        self.schema = dict([("id", "K"), ("ra", "D"), ("dec", "D"), ("thing_id", "K")] +
                           [(f.lower(), "E") for f in FILTERS] +
                           [(f.lower() + "_err", "E") for f in FILTERS])
        self.buildArgs = "-S r -L 20 -E -M -j 0.4 -n 100 -r 1"

    def filter(self, data):
        isGood = data.field("STARNOTGAL").astype("bool")
        return dict([(col, data.field(col.upper())[isGood]) for col in self.schema])
