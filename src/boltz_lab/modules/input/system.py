"""System configuration management for Boltz predictions.

This module provides the System class for managing molecular system
configurations, including proteins, ligands, and their properties.
"""

import yaml
import logging
from boltz_lab.modules.utils import helpers

logger  = logging.getLogger('boltz-lab.helpers')

class System:
    """Manage molecular system configuration for Boltz predictions.

    The System class provides methods to load, update, query, and save
    molecular system definitions from dictionaries or YAML files. System
    configurations typically define proteins, ligands, MSAs, and other
    molecular components.

    Parameters
    ----------
    system : dict, optional
        Pre-loaded system dictionary containing molecular definitions.
    system_path : str, optional
        Path to YAML file to load the system configuration from.

    Raises
    ------
    ValueError
        If both system and system_path are provided, or if neither is provided.

    Examples
    --------
    >>> # Load from dictionary
    >>> system_dict = {"sequences": [{"protein": {"id": "A"}}]}
    >>> sys = System(system=system_dict)

    >>> # Load from YAML file
    >>> sys = System(system_path="system.yaml")
    """
    def __init__(self, system=None, system_path=None):
        self.logger = logging.getLogger('boltz-lab.helpers.system.System')

        if system and system_path:
            raise ValueError("Provide either 'system' or 'system_path', not both.")

        if system:
            self.system = system
        elif system_path:
            self.logger.debug(f"Loading system YAML from {system_path}")
            self.system = helpers.read_yaml(path=system_path)
        else:
            raise ValueError("Either 'system' or 'system_path' must be provided.")
    
    def update_system(self, value, path=None, parent_key=None, sub_key=None):
        """Update system configuration at a specific path or key.

        Modifies the system dictionary either by navigating a specific path
        through nested structures, or by searching for a parent key and
        optionally a sub-key within it.

        Parameters
        ----------
        value : any
            The value to set at the specified location.
        path : list, optional
            Path to the target location as a list of keys/indices.
            Example: ["sequences", 0, "protein", "fasta"]
        parent_key : str, optional
            Top-level or nested key to search for and update.
        sub_key : str, optional
            Sub-key within the parent_key dictionary to update.

        Raises
        ------
        ValueError
            If the path traversal fails due to type mismatch.
        KeyError
            If the parent_key is not found in the system.

        Notes
        -----
        Either path or parent_key must be provided, but not both.
        The path parameter creates nested structures if they don't exist.
        """
        if path:
            d = self.system
            for key in path[:-1]:
                if isinstance(d, dict):
                    d = d.setdefault(key, {})
                elif isinstance(d, list):
                    idx = int(key)
                    d = d[idx]
                else:
                    raise ValueError(f"Cannot traverse into object at {key} in path {path}")

            last_key = path[-1]
            if isinstance(d, dict):
                d[last_key] = value
            elif isinstance(d, list):
                d[int(last_key)] = value
            else:
                raise ValueError(f"Cannot set value at path {path}")

        elif parent_key:
            def set_key(d):
                if isinstance(d, dict):
                    if parent_key in d:
                        if sub_key:
                            d[parent_key] = d.get(parent_key, {})
                            d[parent_key][sub_key] = value
                        else:
                            d[parent_key] = value
                        return True
                    return any(set_key(v) for v in d.values())
                if isinstance(d, list):
                    return any(set_key(i) for i in d)
                return False

            if not set_key(self.system):
                raise KeyError(f"Parent key '{parent_key}' not found in the system.")

    def delete_system_key(self, path=None, parent_key=None, sub_key=None, keys_to_delete=None):
        """
        Delete specific keys from the system configuration at a given path or key.

        Parameters
        ----------
        path : list, optional
            Path to the target dictionary/list as a list of keys/indices.
            Example: ["sequences", 0, "ligand"]
        parent_key : str, optional
            Top-level or nested key to search for and delete keys in.
        sub_key : str, optional
            Sub-key within the parent_key dictionary to target.
        keys_to_delete : list of str, optional
            List of keys to delete from the target dictionary.

        Raises
        ------
        ValueError
            If path traversal fails due to type mismatch.
        KeyError
            If parent_key is not found in the system.
        """
        if keys_to_delete is None:
            return

        if path:
            d = self.system
            for key in path:
                if isinstance(d, dict):
                    d = d.setdefault(key, {})
                elif isinstance(d, list):
                    d = d[int(key)]
                else:
                    raise ValueError(f"Cannot traverse into object at {key} in path {path}")
            if isinstance(d, dict):
                for k in keys_to_delete:
                    if k in d:
                        del d[k]

        elif parent_key:
            def remove_keys(d):
                if isinstance(d, dict):
                    if parent_key in d:
                        target = d[parent_key]
                        if sub_key:
                            target = target.get(sub_key, {})
                            if isinstance(target, dict):
                                for k in keys_to_delete:
                                    target.pop(k, None)
                        else:
                            if isinstance(target, dict):
                                for k in keys_to_delete:
                                    target.pop(k, None)
                        return True
                    return any(remove_keys(v) for v in d.values())
                if isinstance(d, list):
                    return any(remove_keys(i) for i in d)
                return False

            if not remove_keys(sys_obj.system):
                raise KeyError(f"Parent key '{parent_key}' not found in the system.")

    def find_value(self, key=None, path=None):
        """Retrieve a value by path or by recursively searching for a key.

        Supports two modes of retrieval: direct path navigation or recursive
        key search throughout the nested structure.

        Parameters
        ----------
        key : str, optional
            Key name to search for recursively throughout the system.
        path : list, optional
            Specific path to the value as a list of keys/indices.

        Returns
        -------
        any or None
            The value found at the specified location, or None if not found.

        Raises
        ------
        ValueError
            If path navigation fails or if the key appears multiple times.

        Examples
        --------
        >>> sys = System(system={"sequences": [{"protein": {"id": "A", "fasta": "MKRAAT"}}]})

        >>> # Find by path
        >>> sys.find_value(path=["sequences", 0, "protein", "fasta"])
        'MKRAAT'

        >>> # Find by key
        >>> sys.find_value(key="fasta")
        'MKRAAT'

        Notes
        -----
        When using key search, if the key appears multiple times in the
        system, a ValueError is raised. Use path instead for disambiguation.
        """
        if path:
            cur = self.system
            for p in path:
                if isinstance(cur, dict):
                    cur = cur[p]
                elif isinstance(cur, list):
                    cur = cur[int(p)]
                else:
                    raise ValueError(f"Path {path} is invalid")
            return cur

        if key:
            def search(d):
                results = []
                if isinstance(d, dict):
                    for k, v in d.items():
                        if k == key:
                            results.append(v)
                        results.extend(search(v))
                elif isinstance(d, list):
                    for item in d:
                        results.extend(search(item))
                return results

            found = search(self.system)
            if not found:
                logger.info(f"Key not found in system: {key}")
                return None
            if len(found) > 1:
                raise ValueError(f"Key '{key}' appears multiple times; use path instead.")
            return found[0]

    def save_system_to_yaml(self, path):
        """Save the system configuration to a YAML file.

        Writes the current system dictionary to a YAML file, preserving
        the order of keys.

        Parameters
        ----------
        path : str
            Output file path for the YAML file.

        Notes
        -----
        The YAML is written with sort_keys=False to preserve insertion order.
        """
        with open(path, "w") as file:
            yaml.dump(self.system, file, sort_keys=False, default_flow_style=False)