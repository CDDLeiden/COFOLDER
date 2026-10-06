import pandas as pd
import logging

logger = logging.getLogger(__name__)

def parse_censored_affinity(
    affinity_series: pd.Series,
    keep_sign: bool = True
) -> pd.DataFrame:
    """Parse affinity values with censoring signs into numeric and sign components.

    Separates censoring indicators (e.g., '>', '<', '>=', '<=') from the
    numeric affinity values, allowing for downstream analysis of censored data.

    Parameters
    ----------
    affinity_series : pd.Series
        Series of affinity values, possibly as strings with censoring signs
        (e.g., '>5.0', '<=3.2').
    keep_sign : bool, default=True
        Whether to retain the censoring sign in the output DataFrame.
        If False, the 'affinity_sign' column will be set to None.

    Returns
    -------
    pd.DataFrame
        DataFrame with two columns:
        - 'affinity_value' : float - The numeric part of the affinity
        - 'affinity_sign' : str or None - The censoring sign if present

    Examples
    --------
    >>> import pandas as pd
    >>> data = pd.Series(['>5.0', '3.2', '<=2.5'])
    >>> parse_censored_affinity(data)
       affinity_value affinity_sign
    0             5.0             >
    1             3.2          None
    2             2.5            <=
    """
    signs = ['>=', '<=', '>', '<']
    def split_sign(val):
        if pd.isnull(val):
            return (None, None)
        val = str(val).strip()
        for s in signs:
            if val.startswith(s):
                try:
                    return (float(val[len(s):].strip()), s)
                except ValueError:
                    return (None, s)
        try:
            return (float(val), None)
        except ValueError:
            return (None, None)
    parsed = affinity_series.apply(split_sign).tolist()
    df = pd.DataFrame(
        {
            'affinity_value': [value for value, _sign in parsed],
            # Pandas 3 infers a dedicated string dtype here and converts None to
            # NaN.  The public helper promises str-or-None values, so retain an
            # object column explicitly across supported pandas versions.
            'affinity_sign': pd.Series(
                [sign for _value, sign in parsed],
                dtype=object,
            ),
        }
    )
    if not keep_sign:
        df['affinity_sign'] = None
    return df

def remove_censored_affinity(
    df: pd.DataFrame,
    cols: list
) -> pd.DataFrame:
    """Remove rows containing censoring signs in specified columns.

    Filters out rows where any of the specified columns contain censoring
    indicators (>, <, >=, <=), retaining only rows with purely numeric values.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to filter.
    cols : list of str
        List of column names to check for censoring signs.

    Returns
    -------
    pd.DataFrame
        Filtered DataFrame containing only rows where all specified columns
        have numeric values without censoring signs.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'exp': ['5.0', '>3.0', '2.5'],
    ...     'pred': ['4.8', '3.2', '2.3']
    ... })
    >>> remove_censored_affinity(df, ['exp'])
          exp pred
    0     5.0  4.8
    2     2.5  2.3
    """
    import re
    censor_pattern = re.compile(r'^(>=|<=|>|<)')
    mask = pd.Series([True] * len(df))
    for col in cols:
        mask &= ~df[col].astype(str).str.strip().str.match(censor_pattern)
    return df[mask].copy()

def strip_censoring_signs(
    df: pd.DataFrame,
    cols: list
) -> pd.DataFrame:
    """Remove censoring signs from values and convert to floats.

    Strips censoring indicators (>, <, >=, <=) from the beginning of values
    in specified columns and converts them to numeric floats. Non-numeric
    values after stripping are set to NaN.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        List of column names to strip censoring signs from.

    Returns
    -------
    pd.DataFrame
        DataFrame with censoring signs removed and values converted to float
        in the specified columns.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({
    ...     'exp': ['>5.0', '3.0', '<=2.5'],
    ...     'pred': ['4.8', '3.2', '2.3']
    ... })
    >>> strip_censoring_signs(df, ['exp'])
          exp pred
    0     5.0  4.8
    1     3.0  3.2
    2     2.5  2.3
    """
    import re
    df = df.copy()
    censor_pattern = re.compile(r'^(>=|<=|>|<)')
    for col in cols:
        df[col] = df[col].astype(str).str.strip().str.replace(censor_pattern, '', regex=True)
        df[col] = pd.to_numeric(df[col], errors='coerce')
    return df

def prepare_affinity_dataframe(
    df: pd.DataFrame,
    cols: list,
    censoring: str = 'remove'
) -> pd.DataFrame:
    """Prepare a DataFrame for affinity analysis by handling censoring signs.

    Processes a DataFrame to handle censored affinity data, either by
    removing censored rows entirely or by stripping the censoring signs
    and retaining the numeric values.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame to process.
    cols : list of str
        List of column names to check and clean for censoring signs.
    censoring : {'remove', 'strip'}, default='remove'
        Strategy for handling censoring signs:
        - 'remove': Remove rows with censoring signs in any specified column
        - 'strip': Remove censoring signs and use the numeric part

    Returns
    -------
    pd.DataFrame
        Cleaned DataFrame ready for numeric analysis.

    Raises
    ------
    ValueError
        If censoring parameter is not 'remove' or 'strip'.

    Examples
    --------
    >>> import pandas as pd
    >>> df = pd.DataFrame({'exp': ['>5.0', '3.0'], 'pred': ['4.8', '3.2']})
    >>> prepare_affinity_dataframe(df, ['exp'], censoring='remove')
          exp pred
    1     3.0  3.2

    >>> prepare_affinity_dataframe(df, ['exp'], censoring='strip')
          exp pred
    0     5.0  4.8
    1     3.0  3.2
    """
    if censoring == 'remove':
        return remove_censored_affinity(df, cols)
    elif censoring == 'strip':
        return strip_censoring_signs(df, cols)
    else:
        raise ValueError("censoring must be 'remove' or 'strip'")
