"""
Streamlit UI for COFOLDER - A professional interface for running Boltz co-folding workflows.
"""
import os
import sys
import logging
import yaml
from pathlib import Path
from typing import Dict, Any, Optional
import streamlit as st
from streamlit_option_menu import option_menu
import subprocess
from datetime import datetime

from cofolder import __version__

# Add src to path for imports
sys.path.insert(0, str(Path(__file__).parent.parent.parent))

# Configure Streamlit page
st.set_page_config(
    page_title="COFOLDER UI",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# Configure logger
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger("cofolder-ui")

# Default Boltz options
DEFAULT_OPTIONS = {
    'out_dir': './',
    'cache': '~/.boltz',
    'checkpoint': None,
    'devices': 1,
    'accelerator': 'gpu',
    'recycling_steps': 3,
    'sampling_steps': 200,
    'diffusion_samples': 1,
    'max_parallel_samples': 5,
    'step_scale': 1.638,
    'write_full_pae': False,
    'write_full_pde': False,
    'output_format': 'mmcif',
    'num_workers': 2,
    'override': False,
    'seed': None,
    'use_msa_server': False,
    'msa_server_url': 'https://api.colabfold.com',
    'msa_pairing_strategy': 'greedy',
    'use_potentials': False,
    'model': 'boltz2',
    'method': None,
    'affinity_mw_correction': False,
    'sampling_steps_affinity': 200,
    'diffusion_samples_affinity': 5,
    'affinity_checkpoint': None,
    'max_msa_seqs': 8192,
    'subsample_msa': False,
    'num_subsamples_msa': 1024,
    'no_kernels': False,
}

OPTION_GROUPS = {
    'Output & System': [
        'out_dir', 'cache', 'checkpoint', 'override'
    ],
    'Hardware & Performance': [
        'devices', 'accelerator', 'num_workers', 'no_kernels'
    ],
    'Folding Parameters': [
        'recycling_steps', 'sampling_steps', 'diffusion_samples', 'step_scale',
        'max_parallel_samples', 'seed'
    ],
    'Output Formats': [
        'write_full_pae', 'write_full_pde', 'output_format'
    ],
    'MSA Settings': [
        'use_msa_server', 'msa_server_url', 'msa_pairing_strategy',
        'max_msa_seqs', 'subsample_msa', 'num_subsamples_msa'
    ],
    'Potentials & Affinity': [
        'use_potentials', 'affinity_mw_correction', 'sampling_steps_affinity',
        'diffusion_samples_affinity', 'affinity_checkpoint'
    ],
    'Model Settings': [
        'model', 'method'
    ]
}

OPTION_INFO = {
    'out_dir': 'Output directory for results',
    'cache': 'Boltz cache directory for models',
    'checkpoint': 'Path to custom checkpoint file',
    'devices': 'Number of GPU devices to use',
    'accelerator': 'Hardware accelerator (gpu/cpu)',
    'recycling_steps': 'Number of recycling steps in Boltz',
    'sampling_steps': 'Number of sampling steps',
    'diffusion_samples': 'Number of diffusion samples to generate',
    'max_parallel_samples': 'Maximum parallel samples',
    'step_scale': 'Step scaling factor',
    'write_full_pae': 'Write full PAE (Predicted Aligned Error) matrix',
    'write_full_pde': 'Write full PDE matrix',
    'output_format': 'Output structure format (mmcif/pdb)',
    'num_workers': 'Number of worker threads',
    'override': 'Override existing results',
    'seed': 'Random seed for reproducibility',
    'use_msa_server': 'Use MSA server for sequence alignment',
    'msa_server_url': 'URL of MSA server',
    'msa_pairing_strategy': 'Strategy for pairing MSAs',
    'use_potentials': 'Use potential functions',
    'model': 'Model to use (boltz1/boltz2)',
    'method': 'Prediction method',
    'affinity_mw_correction': 'Apply molecular weight correction for affinity',
    'sampling_steps_affinity': 'Sampling steps for affinity prediction',
    'diffusion_samples_affinity': 'Diffusion samples for affinity',
    'affinity_checkpoint': 'Custom checkpoint for affinity model',
    'max_msa_seqs': 'Maximum MSA sequences',
    'subsample_msa': 'Subsample MSA sequences',
    'num_subsamples_msa': 'Number of MSA subsamples',
    'no_kernels': 'Disable kernel optimizations',
}

# Custom CSS for professional styling (theme-aware)
st.markdown("""
    <style>
    :root {
        /* Use Streamlit theme variables only; avoid hard fallbacks so theme switching works */
        --app-bg: var(--background-color);
        --app-bg-2: var(--secondary-background-color);
        --app-text: var(--text-color);
        --app-primary: var(--primary-color);
        --app-accent-2: #764ba2; /* header gradient second color */
        --success-border: rgba(40, 167, 69, 0.6);
        --error-border: rgba(220, 53, 69, 0.7);
        --info-border: rgba(23, 162, 184, 0.6);
    }

    /* Do not force background/text here; let Streamlit theme apply globally */

    .main-header {
        font-size: 2.5em;
        font-weight: bold;
        background: linear-gradient(135deg, var(--app-primary) 0%, var(--app-accent-2) 100%);
        -webkit-background-clip: text;
        -webkit-text-fill-color: transparent;
        background-clip: text;
        margin-bottom: 0.5em;
    }
    .recipe-description {
        font-size: 1.1em;
        color: var(--app-text);
        margin-bottom: 1.5em;
        padding: 1em;
        background-color: var(--app-bg-2);
        border-left: 4px solid var(--app-primary);
        border-radius: 4px;
    }
    .section-header {
        font-size: 1.3em;
        font-weight: bold;
        margin-top: 1.5em;
        margin-bottom: 0.5em;
        padding-bottom: 0.5em;
        border-bottom: 2px solid var(--app-primary);
    }
    .option-container {
        background-color: var(--app-bg-2);
        padding: 1em;
        border-radius: 8px;
        margin-bottom: 0.5em;
        border-left: 3px solid var(--app-primary);
    }
    .success-box, .error-box, .info-box {
        background-color: transparent; /* let theme control background */
        border-radius: 6px;
        padding: 0.9em 1em;
        margin-top: 1em;
        color: var(--app-text);
    }
    .success-box { border: 1px solid var(--success-border); }
    .error-box { border: 1px solid var(--error-border); }
    .info-box { border: 1px solid var(--info-border); }
    </style>
    """, unsafe_allow_html=True)


def initialize_session_state():
    """Initialize Streamlit session state."""
    if 'recipe' not in st.session_state:
        st.session_state.recipe = 'validate'
    if 'options' not in st.session_state:
        st.session_state.options = DEFAULT_OPTIONS.copy()
    if 'execution_log' not in st.session_state:
        st.session_state.execution_log = []


def render_header():
    """Render the main header."""
    col1, col2 = st.columns([3, 1])
    with col1:
        st.markdown('<div class="main-header">🧬 COFOLDER UI</div>', unsafe_allow_html=True)
        st.caption("Professional interface for Boltz protein co-folding workflows")
    with col2:
        st.markdown(f"""
        <div style='text-align: right; margin-top: 1em;'>
            <small>Version: {__version__}</small><br>
            <small>Last updated: {datetime.now().strftime('%Y-%m-%d')}</small>
        </div>
        """, unsafe_allow_html=True)


def load_yaml_file(file_path: str) -> Optional[Dict[str, Any]]:
    """Load and parse a YAML file."""
    try:
        with open(file_path, 'r') as f:
            return yaml.safe_load(f)
    except FileNotFoundError:
        st.error(f"File not found: {file_path}")
        return None
    except yaml.YAMLError as e:
        st.error(f"Error parsing YAML: {e}")
        return None


def save_config(config: Dict[str, Any], file_path: str):
    """Save configuration to a YAML file."""
    try:
        os.makedirs(os.path.dirname(file_path) or '.', exist_ok=True)
        with open(file_path, 'w') as f:
            yaml.dump(config, f, default_flow_style=False, sort_keys=False)
        return True
    except Exception as e:
        st.error(f"Error saving configuration: {e}")
        return False


def render_recipe_selector():
    """Render recipe selection sidebar."""
    st.sidebar.markdown("### 🎯 Select Recipe")

    recipes = {
        'validate': '✅ Validate - Co-fold and validate a single system',
        'screen': '🔬 Screen - Virtual screening library',
        'oracle': '🎰 Oracle - Single SMILES prediction',
    }

    selected = st.sidebar.radio(
        "Choose a recipe:",
        list(recipes.keys()),
        format_func=lambda x: recipes[x],
        key='recipe_selector'
    )

    if selected != st.session_state.recipe:
        st.session_state.recipe = selected
        st.rerun()


def render_validate_ui():
    """Render the Validate workflow UI."""
    st.markdown('<div class="recipe-description">'
                '<b>Validate:</b> Co-fold and validate a single system using Boltz.'
                '</div>', unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown('#### System Configuration')
        system_path = st.text_input(
            'System YAML Path',
            value='./system.yaml',
            help='Path to the system YAML file defining the protein-ligand complex'
        )

    with col2:
        st.markdown('#### Options Configuration')
        options_path = st.text_input(
            'Options YAML Path',
            value='./options.yaml',
            help='Path to the Boltz options YAML file'
        )

    col1, col2 = st.columns(2)
    with col1:
        wrk_dir = st.text_input(
            'Working Directory',
            value=os.getcwd(),
            help='Directory where results will be saved'
        )

    with col2:
        debug = st.checkbox('Enable Debug Logging', value=False)

    # Render options
    render_options_editor(system_path, options_path)

    # Run button
    if st.button('▶️ Run Validation', use_container_width=True, type='primary'):
        run_validate(wrk_dir, system_path, options_path, debug)


def render_screen_ui():
    """Render the Screen (virtual screening) workflow UI."""
    st.markdown('<div class="recipe-description">'
                '<b>Screen:</b> Co-fold a library of compounds for virtual screening. Perfect for large-scale ligand evaluation.'
                '</div>', unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown('#### System Configuration')
        system_path = st.text_input(
            'System YAML Path',
            value='./system.yaml',
            key='screen_system_path',
            help='Path to the system YAML file'
        )

    with col2:
        st.markdown('#### Options Configuration')
        options_path = st.text_input(
            'Options YAML Path',
            value='./options.yaml',
            key='screen_options_path',
            help='Path to the Boltz options YAML file'
        )

    # Ligand replacement configuration
    st.markdown('#### Ligand Replacement')
    col1, col2 = st.columns(2)

    with col1:
        ligand_chain = st.text_input(
            'Ligand Chain',
            value='B',
            help='Existing ligand chain whose chemistry will be replaced'
        )

    with col2:
        library_format = st.selectbox(
            'Library Format',
            ['csv', 'sdf', 'mol'],
            help='CSV uses SMILES columns; SDF/MOL use structure records.'
        )

    # Input file configuration
    col1, col2, col3 = st.columns(3)

    with col1:
        library = st.text_input(
            'Library File Path',
            value='./ligands.csv',
            help='Path to the CSV, SDF, or MOL ligand library'
        )
        smiles_column = (
            st.text_input(
                'SMILES Column',
                value='smiles',
                help='Column containing ligand SMILES'
            )
            if library_format == 'csv'
            else None
        )

    with col2:
        id_field = st.text_input(
            'ID Column' if library_format == 'csv' else 'ID Property',
            value='id' if library_format == 'csv' else '_Name',
            help='Column/property name used for original compound IDs'
        )
        duplicate_id_policy = st.selectbox(
            'Duplicate ID Policy',
            ['reject', 'suffix', 'source_index'],
        )

    with col3:
        generate_conformers = st.selectbox(
            'Generate Conformers',
            [None, '2D', '3D', 'sdf'],
            help='Use source coordinates or generate conformers'
        )

    # Additional columns to merge
    merge_data = st.text_input(
        'Additional Columns to Merge',
        value='',
        help='Comma-separated list of columns/properties to include in output'
    )

    col1, col2 = st.columns(2)
    with col1:
        wrk_dir = st.text_input(
            'Working Directory',
            value=os.getcwd(),
            key='screen_wrk_dir',
            help='Directory where results will be saved'
        )

    with col2:
        debug = st.checkbox('Enable Debug Logging', value=False, key='screen_debug')

    # Render options
    render_options_editor(system_path, options_path, key_suffix='_screen')

    # Run button
    if st.button('▶️ Run Screening', use_container_width=True, type='primary'):
        run_screen(
            wrk_dir, system_path, options_path, ligand_chain,
            library, library_format, smiles_column, id_field,
            duplicate_id_policy, generate_conformers, merge_data, debug
        )


def render_oracle_ui():
    """Render the Oracle workflow UI."""
    st.markdown('<div class="recipe-description">'
                '<b>Oracle:</b> Use Boltz as an oracle function for single SMILES predictions. Great for quick affinity predictions.'
                '</div>', unsafe_allow_html=True)

    col1, col2 = st.columns(2)

    with col1:
        st.markdown('#### System Configuration')
        system_path = st.text_input(
            'System YAML Path',
            value='./system.yaml',
            key='oracle_system_path',
            help='Path to the system YAML file'
        )

    with col2:
        st.markdown('#### Options Configuration')
        options_path = st.text_input(
            'Options YAML Path',
            value='./options.yaml',
            key='oracle_options_path',
            help='Path to the Boltz options YAML file'
        )

    col1, col2 = st.columns(2)
    with col1:
        wrk_dir = st.text_input(
            'Working Directory',
            value=os.getcwd(),
            key='oracle_wrk_dir',
            help='Directory where results will be saved'
        )

    with col2:
        debug = st.checkbox('Enable Debug Logging', value=False, key='oracle_debug')

    # Render options
    render_options_editor(system_path, options_path, key_suffix='_oracle')

    # Run button
    if st.button('▶️ Run Oracle', use_container_width=True, type='primary'):
        run_oracle(wrk_dir, system_path, options_path, debug)


def render_options_editor(system_path: str, options_path: str, key_suffix: str = ''):
    """Render the options editor with grouped controls."""
    st.markdown('#### 🔧 Boltz Options')

    # Load existing options if file exists
    if os.path.exists(options_path):
        loaded_opts = load_yaml_file(options_path)
        if loaded_opts and 'options' in loaded_opts:
            # Convert list of dicts to single dict
            if isinstance(loaded_opts['options'], list):
                for opt_dict in loaded_opts['options']:
                    if isinstance(opt_dict, dict):
                        st.session_state.options.update(opt_dict)

    # Create expandable sections for each group
    for group_name, option_keys in OPTION_GROUPS.items():
        with st.expander(f"📋 {group_name}", expanded=False):
            cols = st.columns(2)
            col_idx = 0

            for option_key in option_keys:
                if option_key not in DEFAULT_OPTIONS:
                    continue

                col = cols[col_idx % 2]
                default_val = st.session_state.options.get(option_key, DEFAULT_OPTIONS[option_key])
                info = OPTION_INFO.get(option_key, '')

                with col:
                    # Render appropriate widget based on option type
                    if isinstance(default_val, bool):
                        st.session_state.options[option_key] = st.checkbox(
                            option_key,
                            value=default_val,
                            help=info,
                            key=f'{option_key}{key_suffix}'
                        )
                    elif isinstance(default_val, int):
                        st.session_state.options[option_key] = st.number_input(
                            option_key,
                            value=default_val,
                            step=1,
                            help=info,
                            key=f'{option_key}{key_suffix}'
                        )
                    elif isinstance(default_val, float):
                        st.session_state.options[option_key] = st.number_input(
                            option_key,
                            value=default_val,
                            step=0.01,
                            help=info,
                            key=f'{option_key}{key_suffix}'
                        )
                    elif option_key in ['model', 'accelerator', 'output_format']:
                        # Special dropdowns
                        if option_key == 'model':
                            options = ['boltz1', 'boltz2']
                        elif option_key == 'accelerator':
                            options = ['gpu', 'cpu']
                        elif option_key == 'output_format':
                            options = ['mmcif', 'pdb']
                        else:
                            options = [str(default_val)]

                        st.session_state.options[option_key] = st.selectbox(
                            option_key,
                            options,
                            index=0 if default_val not in options else options.index(default_val),
                            help=info,
                            key=f'{option_key}{key_suffix}'
                        )
                    else:
                        # String/path input
                        st.session_state.options[option_key] = st.text_input(
                            option_key,
                            value=str(default_val) if default_val is not None else '',
                            help=info,
                            key=f'{option_key}{key_suffix}'
                        )

                col_idx += 1

    # Save options to file
    col1, col2 = st.columns([3, 1])
    with col2:
        if st.button('💾 Save Options', use_container_width=True):
            # Convert options to the expected format
            options_to_save = {
                'wrapper': [
                    {'run_dir': './'},
                    {'system': system_path},
                ],
                'options': [st.session_state.options]
            }
            if save_config(options_to_save, options_path):
                st.success(f"Options saved to {options_path}")


def run_validate(wrk_dir: str, system_path: str, options_path: str, debug: bool):
    """Run the validate workflow."""
    try:
        # Validate inputs
        if not system_path or not options_path:
            st.error('System and Options paths are required.')
            return

        # Save current options
        options_to_save = {
            'version': 1,
            'runtime': {},
            'runner': dict(st.session_state.options),
        }
        save_config(options_to_save, options_path)

        st.info('🚀 Starting validation workflow...')

        # Build command
        cmd = [
            'python', '-m', 'cofolder', 'validate',
            '-w', wrk_dir,
            '-s', system_path,
            '-b', options_path,
        ]

        if debug:
            cmd.append('-d')

        # Run command
        with st.spinner('Running Boltz validation...'):
            result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            st.markdown('<div class="success-box">✅ Validation completed successfully!</div>', unsafe_allow_html=True)
            if result.stdout:
                st.code(result.stdout, language='bash')
        else:
            st.markdown('<div class="error-box">❌ Validation failed!</div>', unsafe_allow_html=True)
            st.code(result.stderr, language='bash')

    except Exception as e:
        st.error(f'Error running validation: {str(e)}')


def run_screen(wrk_dir: str, system_path: str, options_path: str, ligand_chain: str,
               library: Optional[str], library_format: str,
               smiles_column: Optional[str], id_field: str,
               duplicate_id_policy: str, generate_conformers: Optional[str],
               merge_data: Optional[str], debug: bool):
    """Run the screen workflow."""
    try:
        if not all([system_path, options_path, ligand_chain, library, id_field]):
            st.error('System, options, ligand chain, library, and ID field are required.')
            return
        if library_format == 'csv' and not smiles_column:
            st.error('A SMILES column is required for CSV libraries.')
            return

        # Validate input files
        if library and not os.path.exists(library):
            st.error(f'Library file not found: {library}')
            return

        if not library:
            st.error('A compound library is required.')
            return

        # Save current options
        options_to_save = {
            'wrapper': [
                {'run_dir': wrk_dir},
                {'system': system_path},
            ],
            'options': [st.session_state.options]
        }
        save_config(options_to_save, options_path)

        st.info('🚀 Starting screening workflow...')

        # Build command
        cmd = [
            'python', '-m', 'cofolder', 'screen',
            '-w', wrk_dir,
            '-s', system_path,
            '-o', options_path,
            '--ligand_chain', ligand_chain,
            '--library_format', library_format,
            '--duplicate_id_policy', duplicate_id_policy,
            '-c', library,
        ]

        if library_format == 'csv':
            cmd.extend(['--smiles_column', smiles_column, '--col_id', id_field])
        else:
            cmd.extend(['--id_property', id_field])

        if generate_conformers:
            cmd.extend(['--conformers', generate_conformers])

        if merge_data:
            cmd.extend(['--merge_data', merge_data])

        if debug:
            cmd.append('-d')

        # Run command
        with st.spinner('Running virtual screening...'):
            result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            st.markdown('<div class="success-box">✅ Screening completed successfully!</div>', unsafe_allow_html=True)
            if result.stdout:
                st.code(result.stdout, language='bash')
        else:
            st.markdown('<div class="error-box">❌ Screening failed!</div>', unsafe_allow_html=True)
            st.code(result.stderr, language='bash')

    except Exception as e:
        st.error(f'Error running screening: {str(e)}')


def run_oracle(wrk_dir: str, system_path: str, options_path: str, debug: bool):
    """Run the oracle workflow."""
    try:
        if not system_path or not options_path:
            st.error('System and Options paths are required.')
            return

        # Save current options
        options_to_save = {
            'wrapper': [
                {'run_dir': wrk_dir},
                {'system': system_path},
            ],
            'options': [st.session_state.options]
        }
        save_config(options_to_save, options_path)

        st.info('🚀 Starting oracle workflow...')

        # Build command
        cmd = [
            'python', '-m', 'cofolder', 'oracle',
            '-w', wrk_dir,
            '-s', system_path,
            '-b', options_path,
        ]

        if debug:
            cmd.append('-d')

        # Run command
        with st.spinner('Running oracle...'):
            result = subprocess.run(cmd, capture_output=True, text=True)

        if result.returncode == 0:
            st.markdown('<div class="success-box">✅ Oracle completed successfully!</div>', unsafe_allow_html=True)
            if result.stdout:
                st.code(result.stdout, language='bash')
        else:
            st.markdown('<div class="error-box">❌ Oracle failed!</div>', unsafe_allow_html=True)
            st.code(result.stderr, language='bash')

    except Exception as e:
        st.error(f'Error running oracle: {str(e)}')


def main():
    """Main application entry point."""
    initialize_session_state()
    render_header()
    render_recipe_selector()

    # Render selected recipe
    if st.session_state.recipe == 'validate':
        render_validate_ui()
    elif st.session_state.recipe == 'screen':
        render_screen_ui()
    elif st.session_state.recipe == 'oracle':
        render_oracle_ui()

    # Render sidebar information
    with st.sidebar:
        st.markdown('---')
        st.markdown('### 📚 Documentation')
        st.markdown("""
        - **Validate**: Single system co-folding and validation
        - **Screen**: Batch virtual screening
        - **Oracle**: Quick affinity predictions
        """)

        st.markdown('---')
        st.markdown('### ⚙️ Quick Actions')
        if st.button('🔄 Reset Options to Defaults'):
            st.session_state.options = DEFAULT_OPTIONS.copy()
            st.rerun()

        st.markdown('---')
        st.markdown('### 🎨 Appearance')
        st.caption(
            "Use the app menu (⋮ → Settings → Theme) to switch between Light/Dark/Auto. "
            "Your choice is saved per browser.")


if __name__ == '__main__':
    main()
