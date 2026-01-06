Fragment ["(auto)", [3], [3]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212628
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.960615, 1.551458, -1.693306, -2.407707, 2.064065, -2.26086, -3.33713, 2.149564,
                           -1.035076],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([2],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpxfrp460u_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [3], [3]] failed in parallel execution: Fragment calculation failed for ["(auto)", [3], [3]]: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212628
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.960615, 1.551458, -1.693306, -2.407707, 2.064065, -2.26086, -3.33713, 2.149564,
                           -1.035076],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([2],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpxfrp460u_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [7], [7]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212631
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [0.656126, 2.601905, 4.050683, 0.830957, 1.865698, 4.660287, -0.119244, 2.287033, 3.582716],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([6],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp6hl5zwtb_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [13], [13]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212619
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [3.346885, -0.739109, 1.455676, 2.495949, -1.146635, 1.162364, 3.954914, -1.284332, 1.019766],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([12],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpx564rftz_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [5], [5]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212607
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.122137, 3.447414, 1.892218, -2.826831, 4.051366, 2.29064, -1.317729, 3.890921, 2.20323],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([4],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpn9iw3ob9_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [14], [14]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212613
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [0.635713, -1.390592, 3.660639, 1.546539, -1.624232, 3.715803, 0.562533, -1.01904, 2.738487],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([13],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp48hid9t9_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [4], [4]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212608
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [1.294444, 1.23841, 0.012121, 1.076204, 1.09, 0.927564, 1.825869, 0.433479, -0.120713],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([3],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpcrsgq62b_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [8], [8]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212604
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-3.779167, 3.893046, -0.392587, -2.954598, 4.248679, 0.000977, -4.321495, 4.008296, 0.39547],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([7],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp_040bvwi_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [2], [2]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212637
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-0.623906, -0.327068, 1.209419, -0.552412, 0.228582, 0.417315, -1.5621, -0.59689, 1.265256],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([1],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpkw23tqox_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [1], [1]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212620
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [1.289545, 3.622126, -1.234222, 1.516741, 2.736353, -0.934317, 2.006483, 3.786236, -1.881103],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([0],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp0e2c9r_n_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [9], [9]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212616
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [1.517185, -2.080446, 0.275877, 1.370479, -2.993096, -0.031144, 0.767312, -1.634289,
                           -0.177908],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([8],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp419veocn_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [15], [15]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212634
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [4.023571, 2.373326, 2.358002, 4.980223, 2.433122, 2.642029, 3.944691, 1.446061, 2.000753],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([14],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp60f3conu_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [11], [11]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212625
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-3.374181, -1.753947, -1.126278, -3.632256, -1.032976, -1.658281, -3.378912, -2.437496,
                           -1.731673],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([10],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpon9724el_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [12], [12]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212652
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-0.169102, -2.277746, -1.880495, -0.59266, -3.129658, -2.033563, -0.876898, -1.672861,
                           -2.112759],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([11],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmptpcilcw7_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [16], [16]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212649
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [4.15173, -3.364424, 2.461883, 3.843025, -2.655099, 2.963243, 3.347114, -3.828962, 2.183542],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([15],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmppxf98w9u_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [10], [10]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212646
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [1.157642, -3.919991, -4.712879, 1.482782, -4.597069, -5.364193, 1.125355, -3.148726,
                           -5.259923],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([9],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp02ro_7ri_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

Fragment ["(auto)", [6], [6]] execution failed: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212655
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.225668, 2.814883, 4.904774, -2.324615, 1.967662, 5.395445, -1.415078, 3.198525, 5.344975],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([5],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmp6mq1fy96_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

concurrent.futures.process._RemoteTraceback: 
"""
Traceback (most recent call last):
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 321, in _execute_fragment_static
    result = qcng.compute(
             ^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/qcengine/compute.py", line 108, in compute
    output_data = executor.compute(input_data, config)
                  ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/qcengine/programs/psi4.py", line 288, in compute
    raise InputError(error_message)
qcengine.exceptions.InputError: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212628
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.960615, 1.551458, -1.693306, -2.407707, 2.064065, -2.26086, -3.33713, 2.149564,
                           -1.035076],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([2],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpxfrp460u_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.


During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/concurrent/futures/process.py", line 261, in _process_worker
    r = call_item.fn(*call_item.args, **call_item.kwargs)
        ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 712, in _execute_fragment_worker
    label, result = ParallelManyBodyExecutor._execute_fragment_static(task, config)
                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 378, in _execute_fragment_static
    raise RuntimeError(f"Fragment calculation failed for {label}: {exc}")
RuntimeError: Fragment calculation failed for ["(auto)", [3], [3]]: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212628
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.960615, 1.551458, -1.693306, -2.407707, 2.064065, -2.26086, -3.33713, 2.149564,
                           -1.035076],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([2],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpxfrp460u_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

"""

The above exception was the direct cause of the following exception:

Traceback (most recent call last):
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 467, in execute_level_parallel
    result_label, result, duration = future.result()
                                     ^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/concurrent/futures/_base.py", line 449, in result
    return self.__get_result()
           ^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/concurrent/futures/_base.py", line 401, in __get_result
    raise self._exception
RuntimeError: Fragment calculation failed for ["(auto)", [3], [3]]: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212628
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.960615, 1.551458, -1.693306, -2.407707, 2.064065, -2.26086, -3.33713, 2.149564,
                           -1.035076],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([2],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpxfrp460u_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.


During handling of the above exception, another exception occurred:

Traceback (most recent call last):
  File "/home/westh/programming/my_projects/QCManyBody/parallel-execution-project/manuscript_calculations/water_004_qcmanybody.py", line 282, in <module>
    run_parallel_calculation()
  File "/home/westh/programming/my_projects/QCManyBody/parallel-execution-project/manuscript_calculations/water_004_qcmanybody.py", line 236, in run_parallel_calculation
    fragment_results = executor.execute_full_calculation()
                       ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 523, in execute_full_calculation
    level_results = self.execute_level_parallel(level, fragments_at_level)
                    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 473, in execute_level_parallel
    raise RuntimeError(f"Parallel execution failed for fragment {label}: {exc}")
RuntimeError: Parallel execution failed for fragment ["(auto)", [3], [3]]: Fragment calculation failed for ["(auto)", [3], [3]]: 
    -----------------------------------------------------------------------
          Psi4: An Open-Source Ab Initio Electronic Structure Package
                               Psi4 1.10 release

                         Git: Rev {} zzzzzzz 


    D. G. A. Smith, L. A. Burns, A. C. Simmonett, R. M. Parrish,
    M. C. Schieber, R. Galvelis, P. Kraus, H. Kruse, R. Di Remigio,
    A. Alenaizan, A. M. James, S. Lehtola, J. P. Misiewicz, M. Scheurer,
    R. A. Shaw, J. B. Schriber, Y. Xie, Z. L. Glick, D. A. Sirianni,
    J. S. O'Brien, J. M. Waldrop, A. Kumar, E. G. Hohenstein,
    B. P. Pritchard, B. R. Brooks, H. F. Schaefer III, A. Yu. Sokolov,
    K. Patkowski, A. E. DePrince III, U. Bozkaya, R. A. King,
    F. A. Evangelista, J. M. Turney, T. D. Crawford, C. D. Sherrill,
    J. Chem. Phys. 152(18) 184108 (2020). https://doi.org/10.1063/5.0006002

                            Additional Code Authors
    E. T. Seidl, C. L. Janssen, E. F. Valeev, M. L. Leininger,
    J. F. Gonthier, R. M. Richard, H. R. McAlexander, M. Saitow, X. Wang,
    P. Verma, M. H. Lechner, A. Jiang, S. Behnle, A. G. Heide,
    M. F. Herbst, D. L. Poole, T. Győri, E. C. Mitchell, J. P. Pederson,
    and A. M. Wallace

             Previous Authors, Complete List of Code Contributors,
                       and Citations for Specific Modules
    https://github.com/psi4/psi4/blob/master/codemeta.json
    https://github.com/psi4/psi4/graphs/contributors
    https://psicode.org/psi4manual/master/introduction.html#citing-psifour

    -----------------------------------------------------------------------


    Psi4 started on: Tuesday, 30 September 2025 05:01PM

    Process ID: 1212628
    Host:       gamingdesktop
    PSIDATADIR: /home/westh/miniconda3/envs/qcmanybody-parallel/share/psi4
    Memory:     1000.0 MiB
    Threads:    1
    Addons:     a̶d̶c̶c̶, a̶m̶b̶i̶t̶, b̶s̶e̶, c̶c̶t̶3̶, c̶f̶o̶u̶r̶, c̶h̶e̶m̶p̶s̶2̶, c̶p̶p̶e̶, d̶d̶x̶,
                d̶f̶t̶d̶3̶, d̶f̶t̶d̶4̶, dkh, ecpint, einsums, f̶o̶r̶t̶e̶, g̶a̶u̶x̶c̶, g̶c̶p̶,
                g̶d̶m̶a̶, g̶e̶o̶m̶e̶t̶r̶i̶c̶, g̶p̶u̶_̶d̶f̶c̶c̶, i̶n̶t̶e̶g̶r̶a̶t̶o̶r̶x̶x̶,
                i̶p̶i̶, l̶i̶b̶e̶f̶p̶, m̶d̶i̶, m̶p̶2̶d̶, m̶r̶c̶c̶, ooo,
                o̶p̶e̶n̶f̶e̶r̶m̶i̶o̶n̶p̶s̶i̶4̶, pcmsolver, p̶s̶i̶4̶f̶o̶c̶k̶c̶i̶, p̶s̶i̶x̶a̶s̶,
                r̶e̶s̶p̶, s̶i̶m̶i̶n̶t̶, s̶n̶s̶m̶p̶2̶, v̶2̶r̶d̶m̶_̶c̶a̶s̶s̶c̶f̶

  ==> Input QCSchema <==

--------------------------------------------------------------------------
{'driver': 'energy',
 'extras': {},
 'id': None,
 'keywords': {},
 'model': {'basis': 'aug-cc-pvdz', 'method': '(auto)'},
 'molecule': {'extras': {},
              'fix_com': True,
              'fix_orientation': True,
              'fragment_charges': [0.0],
              'fragment_multiplicities': [1],
              'fragments': [[0, 1, 2]],
              'geometry': [-2.960615, 1.551458, -1.693306, -2.407707, 2.064065, -2.26086, -3.33713, 2.149564,
                           -1.035076],
              'masses': [15.99491461957, 1.00782503223, 1.00782503223],
              'molecular_charge': 0.0,
              'molecular_multiplicity': 1,
              'name': 'H32O16 ([2],[])',
              'provenance': {'creator': 'QCElemental',
                             'routine': 'qcelemental.molparse.from_schema',
                             'version': '0.29.0'},
              'real': [True, True, True],
              'schema_name': 'qcschema_molecule',
              'schema_version': 2,
              'symbols': ['O', 'H', 'H'],
              'validated': True},
 'protocols': {},
 'provenance': {'creator': 'QCElemental', 'routine': 'qcelemental.models.results', 'version': '0.29.0'},
 'schema_name': 'qcschema_input',
 'schema_version': 1}
--------------------------------------------------------------------------

Scratch directory: /tmp/tmpxfrp460u_psi_scratch/
Traceback (most recent call last):
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 463, in run_qcschema
    ret_data = run_json_qcschema(input_model.dict(), clean, False, keep_wfn=keep_wfn)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/schema_wrapper.py", line 629, in run_json_qcschema
    val, wfn = methods_dict_[json_data["driver"]](method, **kwargs)
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver.py", line 455, in energy
    plan = task_planner.task_planner("energy", lowername, molecule, **kwargs)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/task_planner.py", line 255, in task_planner
    dermode = negotiate_derivative_type(driver, method, kwargs.pop('dertype', None), verbose=1)
              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 222, in negotiate_derivative_type
    highest_analytic_dertype, proc_messages = highest_analytic_derivative_available(method, proc)
                                              ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/miniconda3/envs/qcmanybody-parallel/lib/python3.11/site-packages/psi4/driver/driver_util.py", line 358, in highest_analytic_derivative_available
    raise MissingMethodError(_alternative_methods_message(method, "any", messages=proc_messages, proc=proc))
psi4.driver.p4util.exceptions.MissingMethodError: Method=(auto) is not available for any derivative level.

