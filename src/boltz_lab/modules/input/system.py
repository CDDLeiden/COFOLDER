import yaml
import logging
from boltz_lab.modules.utils import helpers

logger  = logging.getLogger('boltz-tools.helpers')

class System:
    def __init__(self, system=None, system_path=None):
        """
        Initializes a System object from a dictionary or a YAML file.

        Args:
            system (dict, optional): Pre-loaded system dictionary.
            system_path (str, optional): Path to YAML file to load the system from.
        """
        self.logger = logging.getLogger('boltz-tools.helpers.system.System')

        if system and system_path:
            raise ValueError("Provide either 'system' or 'system_path', not both.")

        if system:
            self.system = system
        elif system_path:
            self.logger.debug(f"Loading system YAML from {system_path}")
            self.system = utils.read_yaml(path=system_path)
        else:
            raise ValueError("Either 'system' or 'system_path' must be provided.")
    
    def update_system(self, value, path=None, parent_key=None, sub_key=None):
        """
        Update system dictionary at a specific path or parent key.

        Args:
            value: Value to set.
            path (list, optional): Path to the value in nested dict/list.
            parent_key (str, optional): Top-level or nested key to update.
            sub_key (str, optional): Sub-key inside parent key.
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

    def find_value(self, key=None, path=None):
        """
        Retrieve a value by path or recursively search by key.

        Args:
            key (str, optional): Key to search recursively.
            path (list, optional): Specific path to the value.
        Returns:
            The value found, or None if not found.
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
        """Save the system dictionary to a YAML file."""
        with open(path, "w") as file:
            yaml.dump(self.system, file, sort_keys=False)