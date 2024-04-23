#ifndef _peephole_hpp_
#define _peephole_hpp_

#include "framework/config.hpp"
#include "framework/utils.hpp"
#include "transpile/circuitopt.hpp"

#include <cassert>
#include <string>

namespace AER {
namespace Transpile {

class Peephole : public CircuitOptimization{
public:
    Peephole() : peephole_enabled_(true) {}
    ~Peephole() {}

    void optimize_circuit(Circuit &circ, Noise::NoiseModel &noise,
                        const opset_t &allowed_opset,
                        ExperimentResult &result) const override;

  // void set_config(const Config &config) override;
  bool enabled() { return peephole_enabled_; }

  //                |0>  |1>   |+>    |->   |L>    |R>    T
  enum BasisState { Low, High, Plus, Minus, Left, Right, Top };
protected:
    mutable bool peephole_enabled_;

    void gate_cancellation(Circuit &circ, const opset_t &allowed_opset) const;
    uint_t check_control_(reg_t &qubits, 
                          std::unordered_map<uint_t, BasisState> &qubitBasisStates,
                          Circuit &circ, uint_t circIdx) const;
    void adjustBasisStates(std::unordered_map<uint_t, BasisState> &qubitBasisState,
                           Operations::Op &op, uint_t opIdx) const;
    void DeleteOperation(Circuit &circ, Operations::Op &op, uint_t opIdx) const;
    void ReplaceOperation(Circuit &circ, uint_t opIdx,
                          bool target, std::string name) const; 

private:
  void dump(const Circuit &circuit) const {
    auto &ops = circuit.ops;
    for (uint_t op_idx = 0; op_idx < ops.size(); ++op_idx) {
      std::cout << std::setw(3) << op_idx << ": ";
      if (ops[op_idx].type == optype_t::nop) {
        std::cout << std::setw(15) << "nop"
                  << ": ";
      } else {
        std::cout << std::setw(15) << ops[op_idx].name << "-"
                  << ops[op_idx].qubits.size() << ": ";
        if (ops[op_idx].qubits.size() > 0) {
          auto qubits = ops[op_idx].qubits;
          std::sort(qubits.begin(), qubits.end());
          int pos = 0;
          for (int j = 0; j < qubits.size(); ++j) {
            int q_pos = 1 + qubits[j] * 2;
            for (int k = 0; k < (q_pos - pos); ++k) {
              std::cout << " ";
            }
            pos = q_pos + 1;
            std::cout << "X";
          }
        }
      }
      std::cout << std::endl;
    }
  }

};

// Transform circ
void Peephole::optimize_circuit(Circuit &circ, Noise::NoiseModel &noise,
                        const opset_t &allowed_opset,
                        ExperimentResult &result) const {
    if (!peephole_enabled_){
        return;
    }
    dump(circ);
    auto &ops = circ.ops;

    gate_cancellation(circ, allowed_opset);
    dump(circ);

}

void Peephole::gate_cancellation(Circuit &circ, const opset_t &allowed_opset) const{
  std::unordered_map<uint_t, BasisState> qubitBasisState;

  // Start all qubits in |0> basis state
  for (uint_t i = 0; i < circ.num_qubits; i++){
    qubitBasisState.insert({i, BasisState::Low});
  }

  for (uint_t i = 0; i < circ.ops.size(); i++){
    // Seperate cx, cz, cy into CNOT and x, y, z
    auto op = circ.ops[i];
    reg_t qubitsOn = circ.ops[i].qubits;

    // 2-Qubit gate operations
    if (qubitsOn.size() > 1){
      // CNOT was removed - Adjust basis as we'll execute the gate
      if (check_control_(qubitsOn, qubitBasisState, circ, i)){
        adjustBasisStates(qubitBasisState, op, i);
      }
    }
    else {
      // 1-Qubit gate operations
      adjustBasisStates(qubitBasisState, op, i);
    }
  }
}

// Return True (1) if the CNOT was removed
uint_t Peephole::check_control_(reg_t &qubits, 
                              std::unordered_map<uint_t, BasisState> &qubitBasisStates,
                              Circuit &circ, uint_t circIdx) const {
  uint_t circuit_changed = 0;
  BasisState controlState = qubitBasisStates[qubits[0]];
  BasisState targetState = qubitBasisStates[qubits[1]];

  switch(controlState){
    case BasisState::Low:
      std::cout << "Erasing at index: " << circIdx << std::endl;
      circ.ops.erase(circ.ops.begin() + circIdx);
      circuit_changed = 1;
      break;
    case BasisState::High:
      if (targetState == BasisState::Top || targetState == BasisState::Low || targetState == BasisState::High){
        // Replace target with X
        ReplaceOperation(circ, circIdx, true, "x");
        circuit_changed = 1;
      } else if (targetState == BasisState::Plus || targetState == BasisState::Minus){
        circ.ops.erase(circ.ops.begin() + circIdx);
        circuit_changed = 1;
      }
      break;
    // Same condition
    case BasisState::Plus:
    case BasisState::Minus:
    case BasisState::Top:
      if (targetState == BasisState::Plus){
        circ.ops.erase(circ.ops.begin() + circIdx);
        circuit_changed = 1;
      } else if (targetState == BasisState::Minus){
        // Replace control with Z
        ReplaceOperation(circ, circIdx, false, "z");
        circuit_changed = 1;
      }
      break;
    default:
      assert(0);
  }
  return circuit_changed;
}

// Basis State Automata
void Peephole::adjustBasisStates(std::unordered_map<uint_t, BasisState> &qubitBasisState,
                                 Operations::Op &op, uint_t opIdx) const {
  if (!op.name.compare("h") && (!op.name.compare("x") || !op.name.compare("cx"))
                           && (!op.name.compare("y") || !op.name.compare("cy")) 
                           && (!op.name.compare("z") || !op.name.compare("cz"))){
    qubitBasisState[opIdx] = BasisState::Top;
    return;
  }
  switch(qubitBasisState[opIdx]){
    case BasisState::Low:
      if (op.name.compare("h")){
        qubitBasisState[opIdx] = BasisState::Plus;
      } else if (op.name.compare("x") || op.name.compare("y")){
        qubitBasisState[opIdx] = BasisState::High;
      }
      break;
    case BasisState::High:
      if (op.name.compare("h")){
        qubitBasisState[opIdx] = BasisState::Minus;
      } else if (op.name.compare("x") || op.name.compare("y")){
        qubitBasisState[opIdx] = BasisState::Low;
      }
      break;
    case BasisState::Plus:
      if (op.name.compare("h")){
        qubitBasisState[opIdx] = BasisState::Low;
      } else if (op.name.compare("z") || op.name.compare("y")){
        qubitBasisState[opIdx] = BasisState::Minus;
      }
      break;
    case BasisState::Minus:
      if (op.name.compare("h")){
        qubitBasisState[opIdx] = BasisState::High;
      } else if (op.name.compare("z") || op.name.compare("y")){
        qubitBasisState[opIdx] = BasisState::Plus;
      }
      break;
    case BasisState::Left:
      if (op.name.compare("h") || op.name.compare("x") || op.name.compare("y")){
        qubitBasisState[opIdx] = BasisState::Right;
      } 
      break;
    case BasisState::Right:
      if (op.name.compare("h") || op.name.compare("x") || op.name.compare("y")){
        qubitBasisState[opIdx] = BasisState::Left;
      } 
      break;
    case BasisState::Top:
      if (op.type == Operations::OpType::reset){
        qubitBasisState[opIdx] = BasisState::Low;
      }
      break;
    default:
      assert(0);
      }
}

void Peephole::ReplaceOperation(Circuit &circ, uint_t opIdx,
                                bool target, std::string name) const {
  Operations::Op rgate;
  rgate.type = Operations::OpType::gate;
  rgate.name = name;
  if (target){
    rgate.qubits = {circ.ops[opIdx].qubits[0]};
  } else{
    rgate.qubits = {circ.ops[opIdx].qubits[1]};
  }
  rgate.string_params = {rgate.name};

  // Insert replacement gate
  circ.ops.emplace(circ.ops.begin() + opIdx, rgate);
  // Remove the CNOT
  circ.ops.erase(circ.ops.begin() + opIdx);
}

} // namespace AER
} // namespace Transpile

#endif
