# KEGG automated metabolite protein interaction network for graph-model (KAMPING)

## Introduction

KEGG features five types of relations: `PPrel`, `GErel`, `PCrel`, `ECrel`, and 
`maplink`. The following figure shows the relation types and their corresponding descriptions.
![img.png](figures/relation_type.png)

Of the five relation types, `ECrel` and `PCrel` describe protein-metabolite interactions. The two entries of `ECrel` 
are two protein (enzyme) entries, with the `value` of the relation being the metabolite entry, it can be `glycan` or 
`compound` (e.g. cpd:C05378 gl:G00037). 

```angular2html
entry1    entry2	type	value	name
hsa:130589	hsa:2538	ECrel	cpd:C00267-90	compound
```



The first entry of `PCrel` is a `compound` entry, and the second entry is a `protein` entry. The `name` and `value` 
of the relation represent the effect of this compound on the protein. The `name` can be `activation`, `inhibition`.

```
entry1    entry2	type	value	name
cpd:C15493-60	hsa:6258	PCrel,PCrel	-->,+p	activation,phosphorylation
```

Due to data parsing, there can be more than one relation between two entries. For example, the following entry has two
the `value` and `name`, the `value` and `name` are separated by a comma.

## Metabolite-protein interaction relation:

We can process `ECrel` relation by expanding it into two binary relation (A-B), also called SIF (simple interaction
format) in BioPAX standard, with first relation with original entry1 as the new entry1 and metabolite as the new
entry2 in the first new relation. Likewise, the second new relation has the original entry2 as the new entry1 and the
metabolite as the new entry2.

```angular2html
entry1    entry2	type	value	name
hsa:130589  cpd:C00267-90 ECrel compound compound
hsa:2538 cpd:C00267-90  ECrel compound compound
# todo: havn't decide the value and name after expanding
```

## Code



After retrieve all relation in an kegg pathway
```
knext mpi --input data/kegg/hsa-ecrel-expanded.txt --output data/kegg/hsa-ecrel-expanded-mpi.txt
```


![img.png](figures/img.png)


