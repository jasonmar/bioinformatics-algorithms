/*
 * Bioinformatics Algorithms in Scala
 * Copyright (C) 2016  Jason Mar
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package bioinf

import java.nio.CharBuffer
import bioinf.Input._
import scala.collection.mutable

object Mutations {
  /**
    * CODE CHALLENGE: Solve the Trie Construction Problem.
    * Input: A collection of strings Patterns.
    * Output: The adjacency list corresponding to Trie(Patterns), in the following format.
    * If Trie(Patterns) has n nodes, first label the root with 0 and then label the remaining nodes with
    * the integers 1 through n - 1 in any order you like. Each edge of the adjacency list of
    * Trie(Patterns) will be encoded by a triple: the first two members of the triple must be the
    * integers labeling the initial and terminal nodes of the edge, respectively; the third member
    * of the triple must be the symbol labeling the edge.

    * Sample Input:
    * ATAGA
    * ATC
    * GAT
    * Sample Output:
    * 0->1:A
    * 1->2:T
    * 2->3:A
    * 3->4:G
    * 4->5:A
    * 2->6:C
    * 0->7:G
    * 7->8:A
    * 8->9:T
    */
  def trieConstruction(patterns: IndexedSeq[String]): Trie = {
    val len = patterns.length * patterns.map{_.length}.max

    val node0 = Array.fill[Int](len){-1}
    val node1 = Array.fill[Int](len){-1}
    val label = Array.fill[Int](len){-1}
    val edges = mutable.Map[(Int,Int),Int]() // edge pointers for (nodeId,label) tuples
    var i = 0 // edge pointer
    var j = 0 // node pointer

    patterns.foreach{pattern =>
      var currentNode = 0 // begin at ROOT node
      pattern.map{label2int}.foreach{currentSymbol =>
        val edge = (currentNode,currentSymbol)
        edges.get(edge) match {
          case Some(idx) => // edge exists
            currentNode = node1(idx)
          case _ => // new edge
            j += 1 // increment node pointer
            node0(i) = currentNode // write initial node to adjacency list
            node1(i) = j // write terminal node to adjacency list
            label(i) = currentSymbol // write label to adjacency list
            currentNode = j // continue from node j for next symbol
            edges += ((edge,i)) // add edge to index
            i += 1 // increment edge pointer
        }
      }

    }
    val data = IndexedSeq(node0,node1,label).map(_.slice(0,i).toIndexedSeq)
    Trie(AdjacencyList(data(0),data(1),data(2)))
  }

  case class AdjacencyList(v: IndexedSeq[Int], w: IndexedSeq[Int], label: IndexedSeq[Int])


  /**
    * CODE CHALLENGE: Implement TRIEMATCHING to solve the Multiple Pattern Matching Problem.
    * Input: A string Text and a collection of strings Patterns.
    * Output: All starting positions in Text where a string from Patterns appears as a substring.

    * TRIEMATCHING(Text, Trie)
    *   while Text is nonempty
    *     PREFIXTRIEMATCHING(Text, Trie)
    *     remove first symbol from Text

    * PREFIXTRIEMATCHING(Text, Trie)
    *   symbol ← first letter of Text
    *   v ← root of Trie
    *   while forever
    *     if v is a leaf in Trie
    *       return the pattern spelled by the path from the root to v
    *     else if there is an edge (v, w) in Trie labeled by symbol
    *       symbol ← next letter of Text
    *       v ← w
    *     else
    *       output "no matches found"
    *       return

    * Sample Input:
    *   AATCGGGTTCAATCGGGGT
    *   ATCG
    *   GGGT
    * Sample Output:
    *   1 4 11 15
    */
  def trieMatching(text: String, patterns: IndexedSeq[String]): IndexedSeq[Int] = {
    val buf = CharBuffer.wrap(text).asReadOnlyBuffer()
    val trie = trieConstruction2(patterns)
    val matchIndices = mutable.ArrayBuffer[Int]()

    var symbol = 0

    while (buf.hasRemaining) {
      symbol = label2int(buf.get())
      buf.mark()
      var v = 0
      var continue = true
      while(continue) {
        val edge = trie.edges0.get((v,symbol))
        val leaf = trie.edges0.get((v,-1)) // pattern termination
        if (leaf.nonEmpty) {
          // match is complete
          buf.reset()
          matchIndices += buf.position() - 1
          continue = false
        } else if (edge.nonEmpty && buf.hasRemaining) {
          symbol = label2int(buf.get()) // next letter of text
          v = trie.adjacencyList.w(edge.get) // next node in trie
        } else {
          // no matches found
          buf.reset()
          continue = false
        }
      }
    }
    matchIndices.result().toIndexedSeq
  }

  case class Trie(adjacencyList: AdjacencyList) {
    lazy val leaves = adjacencyList.w.toSet.diff(adjacencyList.v.toSet)
    lazy val edges = adjacencyList.label.indices.map{i => ((adjacencyList.v(i),adjacencyList.w(i),adjacencyList.label(i)),i)}.toMap
    lazy val edges0 = adjacencyList.label.indices.map{i => ((adjacencyList.v(i),adjacencyList.label(i)),i)}.toMap
    lazy val edges1 = adjacencyList.label.indices.map{i => ((adjacencyList.w(i),adjacencyList.label(i)),i)}.toMap
  }

  /** This trie construction method appends '$' to each pattern to indicate
    * This allows successful resolution of patterns which are nested within each other
    * Example in the book is 'pan' and 'pantry'
    * With terminator characters, they become 'pan$' and 'pantry$' such that 'pan' can be found in the trie
    */
  def trieConstruction2(patterns: IndexedSeq[String]): Trie = {
    trieConstruction(patterns.map{_ + "$"})
  }
}
