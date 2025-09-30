Traceback (most recent call last):
  File "/home/westh/programming/my_projects/QCManyBody/parallel-execution-project/manuscript_calculations/water_002_qcmanybody.py", line 268, in <module>
    run_parallel_calculation()
  File "/home/westh/programming/my_projects/QCManyBody/parallel-execution-project/manuscript_calculations/water_002_qcmanybody.py", line 216, in run_parallel_calculation
    executor = ParallelManyBodyExecutor.from_manybodyinput(
               ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 195, in from_manybodyinput
    return cls(computer_model.qcmb_core, config, driver=driver, specifications=specifications)
           ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
  File "/home/westh/programming/my_projects/QCManyBody/qcmanybody/parallel.py", line 149, in __init__
    self._dependency_graph = core.dependency_graph
                             ^^^^^^^^^^^^^^^^^^^^^
AttributeError: 'ManyBodyCore' object has no attribute 'dependency_graph'
