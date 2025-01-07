from collections import defaultdict
from typing import Any, Optional, Iterable, Union, List, Dict, Literal
import torch
from torch import Tensor


def from_networkx(
        G: Any,
        group_node_attrs: Optional[Union[List[str], Literal['all']]] = None,
        group_edge_attrs: Optional[Union[List[str], Literal['all']]] = None,
) -> 'torch_geometric.data.Data':
    r"""Converts a :obj:`networkx.Graph` or :obj:`networkx.DiGraph` to a
    :class:`torch_geometric.data.Data` instance.

    Args:
        G (networkx.Graph or networkx.DiGraph): A networkx graph.
        group_node_attrs (List[str] or "all", optional): The node attributes to
            be concatenated and added to :obj:`data.x`. (default: :obj:`None`)
        group_edge_attrs (List[str] or "all", optional): The edge attributes to
            be concatenated and added to :obj:`data.edge_attr`.
            (default: :obj:`None`)

    .. note::

        All :attr:`group_node_attrs` and :attr:`group_edge_attrs` values must
        be numeric.

    Examples:
        >>> edge_index = torch.tensor([
        ...     [0, 1, 1, 2, 2, 3],
        ...     [1, 0, 2, 1, 3, 2],
        ... ])
        >>> data = Data(edge_index=edge_index, num_nodes=4)
        >>> g = to_networkx(data)
        >>> # A `Data` object is returned
        >>> from_networkx(g)
        Data(edge_index=[2, 6], num_nodes=4)
    """
    import networkx as nx

    from torch_geometric.data import Data

    G = G.to_directed() if not nx.is_directed(G) else G

    mapping = dict(zip(G.nodes(), range(G.number_of_nodes())))
    edge_index = torch.empty((2, G.number_of_edges()), dtype=torch.long)
    for i, (src, dst) in enumerate(G.edges()):
        edge_index[0, i] = mapping[src]
        edge_index[1, i] = mapping[dst]

    data_dict: Dict[str, Any] = defaultdict(list)
    data_dict['edge_index'] = edge_index

    node_attrs: List[str] = []
    if G.number_of_nodes() > 0:
        node_attrs = list(next(iter(G.nodes(data=True)))[-1].keys())

    edge_attrs: List[str] = []
    if G.number_of_edges() > 0:
        edge_attrs = list(next(iter(G.edges(data=True)))[-1].keys())

    if group_node_attrs is not None and not isinstance(group_node_attrs, list):
        group_node_attrs = node_attrs

    if group_edge_attrs is not None and not isinstance(group_edge_attrs, list):
        group_edge_attrs = edge_attrs

    for i, (_, feat_dict) in enumerate(G.nodes(data=True)):
        if set(feat_dict.keys()) != set(node_attrs):
            raise ValueError('Not all nodes contain the same attributes')
        for key, value in feat_dict.items():
            data_dict[str(key)].append(value)

    for i, (_, _, feat_dict) in enumerate(G.edges(data=True)):
        if set(feat_dict.keys()) != set(edge_attrs):
            raise ValueError('Not all edges contain the same attributes')
        for key, value in feat_dict.items():
            key = f'edge_{key}' if key in node_attrs else key
            data_dict[str(key)].append(value)

    for key, value in G.graph.items():
        if key == 'node_default' or key == 'edge_default':
            continue  # Do not load default attributes.
        key = f'graph_{key}' if key in node_attrs else key
        data_dict[str(key)] = value

    for key, value in data_dict.items():
        if isinstance(value, (tuple, list)) and isinstance(value[0], Tensor):
            data_dict[key] = torch.stack(value, dim=0)
        else:
            try:
                data_dict[key] = torch.as_tensor(value)
            except Exception:
                pass

    # data = Data.from_dict(data_dict)
    data = Data()

    if group_node_attrs is not None:
        xs = []
        for key in group_node_attrs:
            x = data_dict[key]
            x = x.view(-1, 1) if x.dim() <= 1 else x
            xs.append(x)
            # del data[key]
        data.x = torch.cat(xs, dim=-1)

    data.edge_index = data_dict['edge_index']
    if group_edge_attrs is not None:
        xs = []
        for key in group_edge_attrs:
            key = f'edge_{key}' if key in node_attrs else key
            x = data_dict[key]
            x = x.view(-1, 1) if x.dim() <= 1 else x
            xs.append(x)
            # del data[key]
        data.edge_attr = torch.cat(xs, dim=-1)

    if data.x is None and data.pos is None:
        data.num_nodes = G.number_of_nodes()

    return data


def from_hetero_networkx(
        G: Any, node_type_attribute: str,
        edge_type_attribute: Optional[str] = None,
        graph_attrs: Optional[Iterable[str]] = None, nodes: Optional[List] = None,
        group_node_attrs: Optional[Union[List[str], all]] = None,
        group_edge_attrs: Optional[Union[List[str], all]] = None,
        keep_other_attrs: Optional[bool] = False
) -> 'torch_geometric.data.HeteroData':
    r"""Converts a :obj:`networkx.Graph` or :obj:`networkx.DiGraph` to a
    :class:`torch_geometric.data.HeteroData` instance.

    Args:
        G (networkx.Graph or networkx.DiGraph): A networkx graph.
        node_type_attribute (str): The attribute containing the type of a
            node. For the resulting structure to be valid, this attribute
            must be set for every node in the graph. Values contained in
            this attribute will be casted as :obj:`string` if possible. If
            not, the function will raise an error.
        edge_type_attribute (str, optional): The attribute containing the
            type of an edge. If set to :obj:`None`, the value :obj:`"to"`
            will be used in the final structure. Otherwise, this attribute
            must be set for every edge in the graph. (default: :obj:`None`)
        graph_attrs (iterable of str, optional): The graph attributes to be
            copied. (default: :obj:`None`)
        nodes (list, optional): The list of nodes whose attributes are to
            be collected. If set to :obj:`None`, all nodes of the graph
            will be included. (default: :obj:`None`)
        group_node_attrs (List[str] or all, optional): The node attributes to
            be concatenated and added to :obj:`data.x`. They must be present
            for all nodes of each type. (default: :obj:`None`)
        group_edge_attrs (List[str] or all, optional): The edge attributes to
            be concatenated and added to :obj:`data.edge_attr`. They must be
            present for all edge of each type. (default: :obj:`None`)

    Example:

        >>> data = from_hetero_networkx(G, node_type_attribute="type",
        ...                    edge_type_attribute="type")
        <torch_geometric.data.HeteroData()>

    :rtype: :class:`torch_geometric.data.HeteroData`
    """
    import networkx as nx
    from torch_geometric.data import HeteroData

    def get_edge_attributes(G, edge_indexes: list,
                            edge_attrs: list = None) -> dict:
        r"""Collects the attributes of a list of graph edges in a dictionary.

        Args:
            G (networkx.Graph or networkx.DiGraph): A networkx graph.
            edge_indexes (list, optional): The list of edge indexes whose
                attributes are to be collected. If set to :obj:`None`, all
                edges of the graph will be included. (default: :obj:`None`)
            edge_attrs (list, optional): The list of expected attributes to
                be found in every edge. If set to :obj:`None`, the first
                edge encountered will set the values for the rest of the
                process. (default: :obj:`None`)

        Raises:
            ValueError: If some of the edges do not share the same list
            of attributes as the rest, an error will be raised.
        """
        data = defaultdict(list)
        edge_to_data = list(G.edges(data=True))

        for edge_index in edge_indexes:
            _, _, feat_dict = edge_to_data[edge_index]
            if edge_attrs is None:
                edge_attrs = feat_dict.keys()
            if set(feat_dict.keys()) != set(edge_attrs):
                raise ValueError('Not all edges contain the same attributes.')
            for key, value in feat_dict.items():
                data[str(key)].append(value)

        return data

    def get_node_attributes(G, nodes: list,
                            expected_node_attrs: list = None) -> dict:
        r"""Collects the attributes of a list of graph nodes in a dictionary.

        Args:
            G (networkx.Graph or networkx.DiGraph): A networkx graph.
            nodes (list, optional): The list of nodes whose attributes are to
                be collected. If set to :obj:`None`, all nodes of the graph
                will be included. (default: :obj:`None`)
            expected_node_attrs (list, optional): The list of expected
                attributes to be found in every node. If set to :obj:`None`,
                the first node encountered will set the values for the rest
                of the process. (default: :obj:`None`)

        Raises:
            ValueError: If some of the nodes do not share the same
            list of attributes as the rest, an error will be raised.
        """

        data = defaultdict(list)

        node_to_data = G.nodes(data=True)

        for node in nodes:
            feat_dict = node_to_data[node]
            if expected_node_attrs is None:
                expected_node_attrs = feat_dict.keys()
            if set(feat_dict.keys()) != set(expected_node_attrs):
                raise ValueError('Not all nodes contain the same attributes.')
            for key, value in feat_dict.items():
                data[str(key)].append(value)

        return data

    G = G.to_directed() if not nx.is_directed(G) else G

    if nodes is not None:
        G = nx.subgraph(G, nodes)

    hetero_data_dict = {}

    node_to_group_id = {}
    node_to_group = {}
    group_to_nodes = defaultdict(list)
    group_to_edges = defaultdict(list)

    for node, node_data in G.nodes(data=True):
        if node_type_attribute not in node_data:
            raise KeyError(f"Given node_type_attribute: {node_type_attribute} \
                missing from node {node}.")
        node_type = str(node_data[node_type_attribute])
        group_to_nodes[node_type].append(node)
        node_to_group_id[node] = len(group_to_nodes[node_type]) - 1
        node_to_group[node] = node_type

    for i, (node_a, node_b, edge_data) in enumerate(G.edges(data=True)):
        if edge_type_attribute is not None:
            if edge_type_attribute not in edge_data:
                raise KeyError(
                    f"Given edge_type_attribute: {edge_type_attribute} \
                    missing from edge {(node_a, node_b)}.")
            node_type_a, edge_type, node_type_b = edge_data[
                edge_type_attribute]
            if node_to_group[node_a] != node_type_a or node_to_group[
                node_b] != node_type_b:
                raise ValueError(f'Edge {node_a}-{node_b} of type\
                         {edge_data[edge_type_attribute]} joins nodes of types\
                         {node_to_group[node_a]} and { node_to_group[node_b]}.'
                                 )
        else:
            edge_type = "to"
        group_to_edges[(node_to_group[node_a], edge_type,
                        node_to_group[node_b])].append(i)

    for group, group_nodes in group_to_nodes.items():
        hetero_data_dict[str(group)] = {
            k: v
            for k, v in get_node_attributes(G, nodes=group_nodes).items()
            if k != node_type_attribute
        }

    for group, group_edges in group_to_edges.items():
        group_name = '__'.join(group)
        hetero_data_dict[group_name] = {
            k: v
            for k, v in get_edge_attributes(G,
                                            edge_indexes=group_edges).items()
            if k != edge_type_attribute
        }
        edge_list = list(G.edges(data=False))
        global_edge_index = [edge_list[edge] for edge in group_edges]
        group_edge_index = [(node_to_group_id[node_a],
                             node_to_group_id[node_b])
                            for node_a, node_b in global_edge_index]
        hetero_data_dict[group_name]["edge_index"] = torch.tensor(
            group_edge_index, dtype=torch.long).t().contiguous().view(2, -1)

    graph_items = G.graph
    if graph_attrs is not None:
        graph_items = {
            k: v
            for k, v in graph_items.items() if k in graph_attrs
        }
    for key, value in graph_items.items():
        hetero_data_dict[str(key)] = value

    for group, group_dict in hetero_data_dict.items():
        if isinstance(group_dict, dict):
            xs = []
            is_edge_group = group in [
                '__'.join(k) for k in group_to_edges.keys()
            ]
            if is_edge_group:
                group_attrs = group_edge_attrs
            else:
                group_attrs = group_node_attrs
            for key, value in group_dict.items():
                if isinstance(value, (tuple, list)) and isinstance(
                        value[0], torch.Tensor):
                    hetero_data_dict[group][key] = torch.stack(value, dim=0)
                else:
                    try:
                        hetero_data_dict[group][key] = torch.tensor(value)
                    except (ValueError, TypeError):
                        pass
                if group_attrs is not None and key in group_attrs:
                    xs.append(hetero_data_dict[group][key])
            if group_attrs is not None:
                if len(group_attrs) != len(xs):
                    raise KeyError(
                        f'Missing required attribute in group: {group}')
                if is_edge_group:
                    hetero_data_dict[group]['edge_attr'] = torch.cat(
                        xs, dim=-1)
                else:
                    hetero_data_dict[group]['x'] = torch.cat(xs, dim=-1)
        else:
            value = group_dict
            if isinstance(value, (tuple, list)) and isinstance(
                    value[0], torch.Tensor):
                hetero_data_dict[group] = torch.stack(value, dim=0)
            else:
                try:
                    hetero_data_dict[group] = torch.tensor(value)
                except (ValueError, TypeError):
                    pass

    # delete group_node_attrs
    for node_group in group_to_nodes.keys():
        if node_group in hetero_data_dict:
            for node_attr in group_node_attrs:
                del hetero_data_dict[node_group][node_attr]

    if not keep_other_attrs:
        # delete not in group_node_attrs
        for node_group in group_to_nodes.keys():
            if node_group in hetero_data_dict:
                if group_node_attrs is not None:
                    group_node_attrs = group_node_attrs + ['x']
                else:
                    group_node_attrs = ['x']
                keys_to_delete = [node_attr for node_attr in hetero_data_dict[node_group].keys() if node_attr not in group_node_attrs]
                for node_attr in keys_to_delete:
                    del hetero_data_dict[node_group][node_attr]

        # delete not in group_edge_attrs
        for edge_group in group_to_edges.keys():
            edge_group_name = '__'.join(edge_group)
            if edge_group_name in hetero_data_dict:
                if group_edge_attrs is not None:
                    group_edge_attrs = group_edge_attrs + ['edge_attr', 'edge_index']
                else:
                    group_edge_attrs = ['edge_attr', 'edge_index']
                keys_to_delete = [edge_attr for edge_attr in hetero_data_dict[edge_group_name].keys() if edge_attr not in group_edge_attrs]
                for edge_attr in keys_to_delete:
                    del hetero_data_dict[edge_group_name][edge_attr]


    # hetero_data_dict_tensor = {}
    # # delete value in hetero_data_dict if it is not a tensor
    # for key, value in hetero_data_dict.items():
    #     if isinstance(value, dict):
    #         hetero_data_dict_tensor[key] = {}
    #         for k, v in value.items():
    #             if isinstance(v, torch.Tensor):
    #                 hetero_data_dict_tensor[key][k] = v
    #     elif isinstance(value, torch.Tensor):
    #         hetero_data_dict_tensor[key] = value

    return HeteroData(**hetero_data_dict)