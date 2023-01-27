# quekno
A benchmark construction algorithm for quantum circuit transformation

This repository is the Python implementation for the QUEKNO benchmark construction algorithm introduced in 
"On constructing benchmark quantum circuits with near-optimal transformation cost" by Sanjiang Li, Xiangzhen Zhou, and Yuan Feng. https://arxiv.org/abs/2301.08932

Our algorithm uses the subgraph isomorphism algorithm VF2. We construct benchmarks for IBM Q Tokyo (20 qubits), Rochester (53 qubits) and Google's Sycamore (53 qubits). For a different quantum device, you may specify its architecture graph in vf2.ag.py.

The generated benchmarks are saved in ./benchmarkX/ in qasm file, where single-qubit gates are all H gates, and two-qubit gates are all CNOT gate and this is okay for qubit mapping purpose. The near-optimal qct-solutions for these benchmarks are generated and saved in ./metaX/, where a oct-solution for an input circuit is defined as the tuple of an initial mapping and a sequence of actions (each of which is a set of SWAPs).  

Running quekno.py, you can generate QUEKNO benchmarks. What you need to do is to provide AG (architecture graph, can be defined in vf2.ag.py), cost_type, cost, subgraph_size, qbg_ratio etc.

The benchmarks used in the paper (genearted on 22.04.22) are the six zip files like 53Q_depth_Rochester.zip. If you want to read the qct-solution and extract benchmark information such as average cx ratio, you may run benchmark_info_0422.py. You can also find the summarised information in the zip files like meta0422-rochester.zip. 

Send email (Sanjiang Li: mrlisj@gmail.com) if you have questions or suggestions!

@misc{https://doi.org/10.48550/arxiv.2301.08932,
  doi = {10.48550/ARXIV.2301.08932},
  
  url = {https://arxiv.org/abs/2301.08932},
  
  author = {Li, Sanjiang and Zhou, Xiangzhen and Feng, Yuan},
  
  keywords = {Quantum Physics (quant-ph), FOS: Physical sciences, FOS: Physical sciences},
  
  title = {On constructing benchmark quantum circuits with known near-optimal transformation cost},
  
  publisher = {arXiv},
  
  year = {2023},
  
  copyright = {Creative Commons Attribution 4.0 International}
}
