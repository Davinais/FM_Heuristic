#include <iostream>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <cassert>
#include <vector>
#include <cmath>
#include <map>
#include "cell.h"
#include "net.h"
#include "partitioner.h"
using namespace std;


void Partitioner::parseInput(fstream& inFile)
{
    string str;
    // Set balance factor
    inFile >> str;
    _bFactor = stod(str);

    bool initial_part = 0;

    // Set up whole circuit
    while (inFile >> str) {
        if (str == "NET") {
            string netName, cellName, tmpCellName = "";
            inFile >> netName;
            int netId = _netNum;
            _netArray.push_back(new Net(netName));
            _netName2Id[netName] = netId;
            while (inFile >> cellName) {
                if (cellName == ";") {
                    tmpCellName = "";
                    break;
                }
                else {
                    // a newly seen cell
                    if (_cellName2Id.count(cellName) == 0) {
                        int cellId = _cellNum;
                        _cellArray.push_back(new Cell(cellName, 0, cellId));
                        _cellName2Id[cellName] = cellId;
                        _cellArray[cellId]->addNet(netId);
                        _cellArray[cellId]->incPinNum();
                        _netArray[netId]->addCell(cellId);
                        ++_cellNum;
                        tmpCellName = cellName;

                        // Set initial part, to make initial cut size smaller,
                        // we try to put all new cells from a net in the same part
                        _cellArray[cellId]->setPart(initial_part);
                        _netArray[netId]->incPartCount(initial_part);
                        ++_partSize[initial_part];
                    }
                    // an existed cell
                    else {
                        if (cellName != tmpCellName) {
                            assert(_cellName2Id.count(cellName) == 1);
                            int cellId = _cellName2Id[cellName];
                            _cellArray[cellId]->addNet(netId);
                            _cellArray[cellId]->incPinNum();
                            _netArray[netId]->addCell(cellId);
                            tmpCellName = cellName;

                            _netArray[netId]->incPartCount(_cellArray[cellId]->getPart());
                        }
                    }
                }
            }
            initial_part = !initial_part;
            ++_netNum;
        }
    }
    return;
}

void Partitioner::partition()
{
    int min_req_size = static_cast<int>(_cellNum * _bFactor);

    for(auto& net : _netArray) {
        if ((net->getPartCount(0) != 0) && (net->getPartCount(1) != 0)) {
            ++_cutSize;
        }
    }

    do {
        ++_iterNum;
        calc_initial_gain();

        while(_maxGainCell != nullptr) {
            Cell *cell = _cellArray[_maxGainCell->getId()];
            update_gain(_maxGainCell);
        }
        cout << "Best Move " << _bestMoveNum << " with Acc gain " << _maxAccGain << endl;
        cout << "=============================================" << endl;

        while(_moveStack.size() > _bestMoveNum) {
            Cell *mv_cell = _cellArray[_moveStack.back()];
            bool from_part = mv_cell->getPart();
            bool to_part = !from_part;

            mv_cell->move();
            --_partSize[from_part];
            ++_partSize[to_part];
            for(auto& net_id : mv_cell->getNetList()) {
                _netArray[net_id]->decPartCount(from_part);
                _netArray[net_id]->incPartCount(to_part);
            }
            _moveStack.pop_back();
        }
        _moveStack.clear();
        _bList[0].clear();
        _bList[1].clear();

        // _maxAccGain at least 0, so no need to worry if cut size would become larger
        _cutSize -= _maxAccGain;
    } while(_maxAccGain > 0);
}

void Partitioner::calc_initial_gain() {
    _unlockNum[0] = _partSize[0];
    _unlockNum[1] = _partSize[1];

    _moveNum = 0;
    _bestMoveNum = 0;
    _accGain = 0;
    _maxAccGain = 0;

    for(auto& cell : _cellArray) {
        cell->unlock();
        cell->setGain(0);

        bool from_part = cell->getPart();
        bool to_part = !from_part;
        for(auto& net_id : cell->getNetList()) {
            if (_netArray[net_id]->getPartCount(from_part) == 1) {
                cell->incGain();
            }
            if (_netArray[net_id]->getPartCount(to_part) == 0) {
                cell->decGain();
            }
        }
        insert_to_blist(cell->getNode());
    }

    update_max_cell();
}

void Partitioner::update_gain(Node *base_node) {
    Cell *base_cell = _cellArray[base_node->getId()];
    bool from_part = base_cell->getPart();
    bool to_part = !from_part;

    // Move the base cell and lock it
    remove_from_blist(base_node);
    base_cell->move();
    --_partSize[from_part];
    ++_partSize[to_part];
    base_cell->lock();
    --_unlockNum[from_part];

    _moveNum++;
    _moveStack.push_back(base_node->getId());

    _accGain += base_cell->getGain();
    if (_accGain > _maxAccGain) {
        _maxAccGain = _accGain;
        _bestMoveNum = _moveNum;
    }

    for(auto& net_id : base_cell->getNetList()) {
        Net *net = _netArray[net_id];
        net->decPartCount(from_part);
        net->incPartCount(to_part);

        int from_count = net->getPartCount(from_part);
        int to_count = net->getPartCount(to_part);

        for(auto& cell_id : net->getCellList()) {
            Cell *cell = _cellArray[cell_id];
            if (cell->getLock() == false) {
                remove_from_blist(cell->getNode());

                if (to_count == (0 + 1)) {
                    cell->incGain();
                }
                else if (to_count == (1 + 1)) {
                    if (cell->getPart() == to_part) {
                        cell->decGain();
                    }
                }

                if (from_count == 0) {
                    cell->decGain();
                }
                else if (from_count == 1) {
                    if (cell->getPart() == from_part) {
                        cell->incGain();
                    }
                }

                insert_to_blist(cell->getNode());
            }
        }
    }

    update_max_cell();
}

void Partitioner::update_max_cell() {
    int min_req_size = static_cast<int>(_cellNum * _bFactor);

    auto it_0 = _bList[0].rbegin();
    auto it_1 = _bList[1].rbegin();
    if ((!_bList[0].empty()) && (!_bList[1].empty())) {
        if (it_0->first >= it_1->first) {
            _maxGainCell = (_partSize[0] > min_req_size) ? it_0->second : it_1->second;
        }
        else {
            _maxGainCell = (_partSize[1] > min_req_size) ? it_1->second : it_0->second;
        }
    }
    else {
        if ((_bList[0].empty()) && (!_bList[1].empty())) {
            _maxGainCell = (_partSize[1] > min_req_size) ? it_1->second : nullptr;
        }
        else if ((!_bList[0].empty()) && (_bList[1].empty())) {
            _maxGainCell = (_partSize[0] > min_req_size) ? it_0->second : nullptr;
        }
        else {
            _maxGainCell = nullptr;
        }
    }
}

void Partitioner::insert_to_blist(Node *ins_node) {
    Cell *ins_cell = _cellArray[ins_node->getId()];
    auto ret = _bList[ins_cell->getPart()].insert({ins_cell->getGain(), ins_node});
    if (ret.second == false) {
        Node *current = ret.first->second;
        if (current->getNext()) {
            current->getNext()->setPrev(ins_node);
        }
        ins_node->setNext(current->getNext());
        ins_node->setPrev(current);
        current->setNext(ins_node);
    }
}

void Partitioner::remove_from_blist(Node *rm_node) {
    Cell *rm_cell = _cellArray[rm_node->getId()];
    if (rm_node->getPrev() == nullptr) {
        if (rm_node->getNext() == nullptr) {
            _bList[rm_cell->getPart()].erase(rm_cell->getGain());
        }
        else {
            _bList[rm_cell->getPart()][rm_cell->getGain()] = rm_node->getNext();
            rm_node->getNext()->setPrev(nullptr);
        }
    }
    Node::unlink(rm_node);
}

void Partitioner::printSummary() const
{
    cout << endl;
    cout << "==================== Summary ====================" << endl;
    cout << " Cutsize: " << _cutSize << endl;
    cout << " Total cell number: " << _cellNum << endl;
    cout << " Total net number:  " << _netNum << endl;
    cout << " Cell Number of partition A: " << _partSize[0] << endl;
    cout << " Cell Number of partition B: " << _partSize[1] << endl;
    cout << "=================================================" << endl;
    cout << endl;
    return;
}

void Partitioner::reportNet() const
{
    cout << "Number of nets: " << _netNum << endl;
    for (size_t i = 0, end_i = _netArray.size(); i < end_i; ++i) {
        cout << setw(8) << _netArray[i]->getName() << ": ";
        vector<int> cellList = _netArray[i]->getCellList();
        for (size_t j = 0, end_j = cellList.size(); j < end_j; ++j) {
            cout << setw(8) << _cellArray[cellList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::reportCell() const
{
    cout << "Number of cells: " << _cellNum << endl;
    for (size_t i = 0, end_i = _cellArray.size(); i < end_i; ++i) {
        cout << setw(8) << _cellArray[i]->getName() << ": ";
        vector<int> netList = _cellArray[i]->getNetList();
        for (size_t j = 0, end_j = netList.size(); j < end_j; ++j) {
            cout << setw(8) << _netArray[netList[j]]->getName() << " ";
        }
        cout << endl;
    }
    return;
}

void Partitioner::writeResult(fstream& outFile)
{
    stringstream buff;
    buff << _cutSize;
    outFile << "Cutsize = " << buff.str() << '\n';
    buff.str("");
    buff << _partSize[0];
    outFile << "G1 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 0) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    buff.str("");
    buff << _partSize[1];
    outFile << "G2 " << buff.str() << '\n';
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        if (_cellArray[i]->getPart() == 1) {
            outFile << _cellArray[i]->getName() << " ";
        }
    }
    outFile << ";\n";
    return;
}

void Partitioner::clear()
{
    for (size_t i = 0, end = _cellArray.size(); i < end; ++i) {
        delete _cellArray[i];
    }
    for (size_t i = 0, end = _netArray.size(); i < end; ++i) {
        delete _netArray[i];
    }
    return;
}
