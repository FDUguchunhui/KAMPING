# KGML types
# for relation, only gene, compound, group, map
# for species-specfic 'enzyme', 'reaction', 'reaction', 'brite', 'other' is not included
ENTRY_TYPES = [ 'ortholog', 'gene', 'group', 'compound', 'map']
# ignore  'maplink' for now
RELATION_TYPES = ['PPrel', 'GErel', 'PCrel', 'ECrel']

COMPOUND_PROPAGATION_TYPES = 'PPrel'
COMPOUND_PROPAGATION_NAME = 'compound-propagation'
COMPOUND_PROPAGATION_VALUE = 'custom'
COMPOUND_PROPAGATION_DIRECTION = 'directed'

GENE_PROPAGATION_TYPES = 'PPrel'
GENE_PROPAGATION_NAME = 'gene-propagation'
GENE_PROPAGATION_VALUE = 'custom'
GENE_PROPAGATION_DIRECTION = 'directed'

GROUP_EXPANSION_TYPES = 'PPrel'
GROUP_EXPANSION_NAME = 'group-expansion'
GROUP_EXPANSION_VALUE = 'custom'
GROUP_EXPANSION_DIRECTION = 'undirected'

MULTI_SUBSTRATE_TYPES = 'CCrel'
MULTI_SUBSTRATE_NAME = 'multi-substrate'
MULTI_SUBSTRATE_VALUE = 'custom'
MULTI_SUBSTRATE_DIRECTION = 'undirected'

MULTI_PRODUCT_TYPES = 'CCrel'
MULTI_PRODUCT_NAME = 'multi-product'
MULTI_PRODUCT_VALUE = 'custom'
MULTI_PRODUCT_DIRECTION = 'undirected'