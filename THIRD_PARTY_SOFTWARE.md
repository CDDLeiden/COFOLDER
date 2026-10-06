# Third-party software

The COFOLDER license applies only to software and other materials for which the
authors hold the necessary rights. It does not relicense the third-party
software below. Users remain responsible for complying with each upstream
license and its terms.

| Resource | Files or role in this repository | Source/version | Identifier | Verified license | Notes |
|---|---|---|---|---|---|
| SuCOS | `src/cofolder/modules/analytics/sucos.py` contains a minimal adapted implementation | Upstream `master`; the exact revision used for the adaptation was not recorded | [upstream repository](https://github.com/susanhleung/SuCOS) | [MIT License](https://github.com/susanhleung/SuCOS/blob/2cc14509199504719d33964aaa0a9e01c90d6f97/LICENSE.txt), verified at commit `2cc14509199504719d33964aaa0a9e01c90d6f97` | COFOLDER implements the SuCOS shape and pharmacophore-feature scoring concept; see the notice below. |
| MMseqs2 | External binary used by `src/cofolder/tools/fetch_bias_training_data.py` and `src/cofolder/tools/build_bias_training_data.py`; optionally installed by `src/cofolder/tools/install_mmseqs.py`. Corresponding checkout scripts delegate to these installed tools | The optional installer downloads the upstream `latest` Linux AVX2 archive and verifies its configured SHA-256; an MMseqs2 release is not otherwise pinned | [upstream repository](https://github.com/soedinglab/MMseqs2) | [MIT License](https://github.com/soedinglab/MMseqs2/blob/eec9c354be4276d2373996af2e50808b1390d527/LICENSE.md), verified at commit `eec9c354be4276d2373996af2e50808b1390d527` | MMseqs2 is invoked as an external program; its source code is not copied into COFOLDER. |
| boltztools | `src/cofolder/modules/entities/ligand.py` (`mol_to_ccd`) contains an adapted implementation | Upstream `main`; the exact revision used for the adaptation was not recorded | [upstream repository](https://github.com/jacktday/boltztools) | [MIT License](https://github.com/jacktday/boltztools/blob/8dd611920f7464bab723756ff4b4efa93e5c692c/LICENSE), verified at commit `8dd611920f7464bab723756ff4b4efa93e5c692c` | The adapted function prepares and caches RDKit molecules in the Boltz CCD format; see the notice below. |

## SuCOS attribution and adaptation notice

`src/cofolder/modules/analytics/sucos.py` contains a minimal implementation
adapted from [susanhleung/SuCOS](https://github.com/susanhleung/SuCOS). It uses
RDKit shape protrusion distance and pharmacophore feature-map overlap to
calculate the composite SuCOS score. The exact upstream revision used when the
implementation was adapted was not recorded. The upstream MIT license was
verified at the immutable commit linked in the table above.

### MIT License

Copyright &lt;2019&gt; &lt;University of Oxford&gt;

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## MMseqs2 attribution and usage notice

COFOLDER invokes [MMseqs2](https://github.com/soedinglab/MMseqs2) as an
external binary when fetching or building protein bias-training data. The
optional vendor installer downloads a precompiled binary but does not add that
binary to this repository. Because the installer uses an upstream `latest`
URL, the software release is not pinned by that URL; the installer instead
checks the downloaded archive against its configured SHA-256. The upstream MIT
license was verified at the immutable commit linked in the table above.

### MIT License

Copyright (c) 2024 The MMseqs2 Development Team

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.

## boltztools attribution and adaptation notice

`src/cofolder/modules/entities/ligand.py` contains a `mol_to_ccd` implementation
adapted from [jacktday/boltztools](https://github.com/jacktday/boltztools). The
function prepares RDKit molecules, derives geometric and stereochemical
constraints, and caches the result for Boltz. The exact upstream revision used
when the implementation was adapted was not recorded. The upstream MIT license
was verified at the immutable commit linked in the table above.

### MIT License

Copyright (c) 2025 Jack Day

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
