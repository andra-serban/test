// Copyright 2019 Andra Serban, Nicolae Mara
#include <vector>
#include <string>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <set>
#include <map>
#include <utility>
#include <iterator>
#include <cassert>
#include "./aegraph.h"

template <typename T>
void removeDuplicates(std::vector<T> &v) {
    auto end = v.end();
    for (auto it = v.begin(); it != end; ++it) {
        end = std::remove(it + 1, end, *it);
    }
    v.erase(end, v.end());
}

std::string strip(std::string s) {
    // deletes whitespace from the beginning and end of the string
    s.erase(0, s.find_first_not_of(" \n\r\t"));
    s.erase(s.find_last_not_of(" \n\r\t")+1);
    return s;
}

std::pair<std::string, std::string> split_first(std::string s,
    char delimiter = ',') {
    // returns a pair that contains: <first_cut, rest_of_graph>

    int numOpen = 0;

    int len = s.size();
    for (int i = 0; i < len; ++i) {
        char c = s[i];
        if (c == delimiter && numOpen == 0)
            return std::make_pair(
                    strip(s.substr(0, i)), strip(s.substr(i+1, len)));
        if (c == '[')
            numOpen += 1;
        if (c == ']')
            numOpen -= 1;
    }

    return {strip(s), std::string()};
}


std::vector<std::string> split_level(std::string s, char delimiter = ',') {
    // splits 's' into separate entities (atoms, subgraphs)

    std::vector<std::string> result;
    auto snd = s;
    while (true) {
        auto p = split_first(snd, delimiter);
        auto fst = p.first;
        snd = p.second;

        result.push_back(fst);

        if (snd.empty())
            return result;
    }
}


int AEGraph::num_subgraphs() const {
    return subgraphs.size();
}


int AEGraph::num_atoms() const {
    return atoms.size();
}


int AEGraph::size() const {
    return num_atoms() + num_subgraphs();
}


bool AEGraph::operator<(const AEGraph& other) const {
    return this->repr() < other.repr();
}

bool AEGraph::operator==(const AEGraph& other) const {
    return this->repr() == other.repr();
}

bool AEGraph::operator!=(const AEGraph& other) const {
    return this->repr() != other.repr();
}

AEGraph AEGraph::operator[](const int index) const {
    // offers an easier way of accessing the nested graphs
    if (index < num_subgraphs()) {
        return subgraphs[index];
    }

    if (index < num_subgraphs() + num_atoms()) {
        return AEGraph('(' + atoms[index - num_subgraphs()] + ')');
    }

    return AEGraph("()");
}

std::ostream& operator<<(std::ostream &out, const AEGraph &g) {
    out << g.repr();
    return out;
}

AEGraph::AEGraph(std::string representation) {
    // constructor that creates an AEGraph structure from a
    // serialized representation
    char left_sep = representation[0];
    char right_sep = representation[representation.size() - 1];

    assert((left_sep == '(' && right_sep == ')')
        || (left_sep == '[' && right_sep == ']'));

    // if the left separator is '(' then the AEGraph is the entire
    // sheet of assertion
    if (left_sep == '(') {
        is_SA = true;
    } else {
        is_SA = false;
    }

    // eliminate the first pair of [] or ()
    representation = representation.substr(1, representation.size() - 2);

    // split the graph into separate elements
    auto v = split_level(representation);
    // add them to the corresponding vector
    for (auto s : v) {
        if (s[0] != '[') {
            atoms.push_back(s);
        } else {
            subgraphs.push_back(AEGraph(s));
        }
    }

    // also internally sort the new graph
    this->sort();
}

std::string AEGraph::repr() const {
    // returns the serialized representation of the AEGraph

    std::string left, right;
    if (is_SA) {
        left = '(';
        right = ')';
    } else {
        left = '[';
        right = ']';
    }

    std::string result;
    for (auto subgraph : subgraphs) {
        result += subgraph.repr() + ", ";
    }

    int len = atoms.size();
    if (len != 0) {
        for (int i = 0; i < len - 1; ++i) {
            result += atoms[i] + ", ";
        }
        result += atoms[len - 1];
    } else {
        if (subgraphs.size() != 0)
            return left + result.substr(0, result.size() - 2) + right;
    }

    return left + result + right;
}


void AEGraph::sort() {
    std::sort(atoms.begin(), atoms.end());

    for (auto& sg : subgraphs) {
        sg.sort();
    }

    std::sort(subgraphs.begin(), subgraphs.end());
}

bool AEGraph::contains(const std::string other) const {
    // checks if an atom is in a graph
    if (find(atoms.begin(), atoms.end(), other) != atoms.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

bool AEGraph::contains(const AEGraph& other) const {
    // checks if a subgraph is in a graph
    if (find(subgraphs.begin(), subgraphs.end(), other) != subgraphs.end())
        return true;

    for (const auto& sg : subgraphs)
        if (sg.contains(other))
            return true;

    return false;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const std::string other)
    const {
    // returns all paths in the tree that lead to an atom like <other>
    std::vector<std::vector<int>> paths;

    int len_atoms = num_atoms();
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_atoms; ++i) {
        if (atoms[i] == other && size() > 1) {
            paths.push_back({i + len_subgraphs});
        }
    }

    for (int i = 0; i < len_subgraphs; ++i) {
        if (subgraphs[i].contains(other)) {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }

    return paths;
}

std::vector<std::vector<int>> AEGraph::get_paths_to(const AEGraph& other)
    const {
    // returns all paths in the tree that lead to a subgraph like <other>
    std::vector<std::vector<int>> paths;
    int len_subgraphs = num_subgraphs();

    for (int i = 0; i < len_subgraphs; ++i) {
        if (subgraphs[i] == other && size() > 1) {
            paths.push_back({i});
        } else {
            auto r = subgraphs[i].get_paths_to(other);
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
        }
    }
    return paths;
}

std::vector<std::vector<int>> AEGraph::possible_double_cuts() const {
    // se verifica daca se pot gasi doua negatii consecutive
    std::vector<std::vector<int>> paths;
    for (int i = 0; i < num_subgraphs(); ++i) {
        if (this->subgraphs[i].num_subgraphs() == 1 &&
            this->subgraphs[i].num_atoms() == 0) {
            // se adauga fiecare possible double cut intr-o linie din matrice
            paths.push_back({i});
        }
    }
    for (int i = 0; i < num_subgraphs(); ++i) {
        // se parcurg recursiv subgrafurile
         auto r = subgraphs[i].possible_double_cuts();
            for (auto& v : r)
                v.insert(v.begin(), i);
            copy(r.begin(), r.end(), back_inserter(paths));
    }
    return paths;
}

AEGraph AEGraph::double_cut(std::vector<int> where) const {
    // noul graf ce va fi returnat
    AEGraph newGraph = *this;
    // pointer prin care se parcurge graful
    AEGraph *traversalNode = &newGraph;
    int i, pathLength = where.size() - 1;
    // se parcurge graful pana la nodul unde se gaseste double cut
    for (i = 0; i < pathLength; ++i) {
        traversalNode = &traversalNode->subgraphs[where[i]];
    }
    // in child se retine nodul caruia i s-au eliminat cele doua negatii
    AEGraph child = traversalNode->subgraphs[where[i]];
    child = child.subgraphs[0];
    // se sterge subgraful vechi ce contine cele doua negatii
    traversalNode->subgraphs.erase(traversalNode->subgraphs.begin() +
        where[i]);
    // atomii si sugrafurile acestui nod se concateneaza vectorilor atoms si
    // subgraphs ale lui traversalNode
    for (unsigned int i = 0 ; i < child.atoms.size(); ++i) {
        traversalNode->atoms.push_back(child.atoms[i]);
    }

    for (unsigned int i = 0 ; i < child.subgraphs.size(); ++i) {
        traversalNode->subgraphs.push_back(child.subgraphs[i]);
    }
    // am incercat cu copy si back_inserter dar nu  mers (・о・)
    return newGraph;
}

std::vector<std::vector<int>> AEGraph::possible_erasures(int level) const {
    std::vector<std::vector<int>> paths;
    // daca copiii se afla pe nivel par si numarul de subgrafuri este 0, dar
    // exista un atom, respectiv daca numarul de atomi este 0, dar exista un
    // subgraf, se adauga in pat
    if (level % 2 != 0 && !(num_subgraphs() == 0 && num_atoms() == 1 && level
     != -1) && !(num_subgraphs() == 1 && num_atoms() == 0 && level != -1)) {
        for (int i = 0; i < num_subgraphs() + num_atoms(); ++i) {
            paths.push_back({i});
        }
    }
    // apelam recursiv functia pentru fiecare subgraf
    for (int i = 0; i < num_subgraphs(); ++i) {
        std::vector<std::vector<int>> r;
        r = subgraphs[i].possible_erasures(level + 1);
        for (auto& v : r) {
                v.insert(v.begin(), i);
        }
        copy(r.begin(), r.end(), back_inserter(paths));
    }
    return paths;
}

AEGraph AEGraph::erase(std::vector<int> where) const {
    // noul graf ce va fi returnat
    AEGraph newGraph = *this;
    // pointer prin care se parcurge graful
    AEGraph *traversalNode = &newGraph;
    int i, pathLength = where.size() - 1;
    // se parcurge graful pana la nodul unde se gaseste pozitia pe care se face
    // erase
    // restul pasilor sunt asemanatori cu cei din functia deiterate
    for (i = 0; i < pathLength; ++i) {
        traversalNode = &traversalNode->subgraphs[where[i]];
    }
    if (traversalNode->num_subgraphs() == 0) {
        traversalNode->atoms.erase(traversalNode->atoms.begin() +
            where[pathLength]);
    } else if (traversalNode->num_subgraphs() - 1 >= where[pathLength]) {
        traversalNode->subgraphs.erase(traversalNode->subgraphs.begin() +
            where[pathLength]);
    } else {
        traversalNode->atoms.erase(traversalNode->atoms.begin() +
            where[pathLength] - num_subgraphs());
    }
    return newGraph;
}

std::vector<std::vector<int>> AEGraph::possible_deiterations() const {
    std::vector<std::vector<int>> paths;
    // se cauta in ceilalti descendenti ai nodului curent atomii nodului
    for (unsigned int i = 0; i < atoms.size(); ++i) {
        for (unsigned int j = 0; j < subgraphs.size(); ++j) {
            auto aux = subgraphs[j].get_paths_to(atoms[i]);
            for (auto &v : aux) {
                v.insert(v.begin(), j);
            }
            copy(aux.begin(), aux.end(), back_inserter(paths));
        }
    }
    // se cauta in ceilalti descendenti ai nodului curent subgraful fiu al
    // nodului
    for (unsigned int i = 0; i < subgraphs.size(); ++i) {
        for (unsigned int j = 0; j < subgraphs.size(); ++j) {
            if (i != j) {
                auto aux = subgraphs[j].get_paths_to(subgraphs[i]);
                for (auto &v : aux) {
                    v.insert(v.begin(), j);
                }
                copy(aux.begin(), aux.end(), back_inserter(paths));
            }
        }

        auto aux = subgraphs[i].possible_deiterations();
        for (auto &v : aux) {
            v.insert(v.begin(), i);
        }
        copy(aux.begin(), aux.end(), back_inserter(paths));
    }
    // se elimina dublurile
    removeDuplicates(paths);

    return paths;
}

AEGraph AEGraph::deiterate(std::vector<int> where) const {
    // noul graf ce va fi returnat
    AEGraph newGraph = *this;
    // pointer prin care se parcurge graful
    AEGraph *traversalNode = &newGraph;
    int i, pathLength = where.size() - 1;
    // se parcurge graful pana la nodul unde se gaseste double cut
    for (i = 0; i < pathLength; ++i) {
        traversalNode = &traversalNode->subgraphs[where[i]];
    }

    // trebuie se stearaga copilul de pe pozitia where[i]
    // dar trebuie verificat daca este subgraf sau atom
    // am observat ca in graf:
    // 1. subgrafurile sunt reprezentate in partea stanga
    // 2. atomii sunt reprezentati in partea dreapta
    if ((int)traversalNode->subgraphs.size() > where[i]) {
        traversalNode->subgraphs.erase(traversalNode->subgraphs.begin()
            + where[i]);
        return newGraph;
    } else {
        traversalNode->atoms.erase(traversalNode->atoms.begin() + (where[i] -
            traversalNode->subgraphs.size()));
    }
    return newGraph;
}

