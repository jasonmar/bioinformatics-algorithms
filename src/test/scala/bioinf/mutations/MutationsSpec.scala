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

package bioinf.mutations

import bioinf.UnitSpec
import bioinf.mutations.Mutations._

class MutationsSpec extends UnitSpec {
  "Mutations" should "locate pattern in SuffixTree" in {
    val text = "CCAAGCTGCTAGAGG"
    val pattern = "AA"
    val s = suffixTreeConstruction(text)
    patternExistsInSuffixTree(pattern,s) should be (right = true)

    val text4 = "CACCTCGTCAATACAACAAAAGGCGGCTCGCTTAAAGGGCGCAGCTAGTTCCTCCCCCTCTCATTGGGACATAGTCAACCTGCTAATCCGGATTCGAATGGATTATTCCGTAATTGAACGGTAATTTAGTGAGCTTCGCAGTAAACGATAGATGCGAGCTCTAGCAGGCCACTGACTATATAAACGCCAACACTAGTGCCGTGCATGGACGACTCGATGTACTATAGATTTGCACAGGTATGACCGGAGGAGCGGGACTGCCTAGGCTATAGGGAACGGGGAGTATTGGGAGCCTTTTAGGCCCTCGTCATATCCCTTAACGTTCCCGCGCAGCTAAATTGTGGAACCGGAAAACAATGGATCTGCTTATTTTTGTAGGCTTGGTTAAGCGAAACGGATCAAAATAAACAAAGAATTAATCAATGAACTAACCAACGAAGTAAGCAAGGATATACATAGATTTATTCATTGATCTATCCATCGATGTATGCATGGACACAGACTTACTCACTGACCTACCCACCGACGTACGCACGGAGAGTTAGTCAGTGAGCTAGCCAGCGAGGTAGGCAGGGTTTTCTTTGTTCCTTCGTTGCTTGGTCTCTGTCCCTCCGTCGCTCGGTGTGCCTGCGTGGCTGGGCCCCGCCGGCGCGGGGAAA"
    val s2 = suffixTreeConstruction(text4)

    patternExistsInSuffixTree("AA",s2) should be (right = true)
    patternExistsInSuffixTree("AAA",s2) should be (right = true)
    patternExistsInSuffixTree("GATT",s2) should be (right = true)
    patternExistsInSuffixTree("GGATT",s2) should be (right = true)
    patternExistsInSuffixTree("AAACC",s2) should be (right = false)
  }

  it should "identify the shortest non-shared substring" in {
    val text1 = "CCAAGCTGCTAGAGG"
    val text2 = "CATGCTGGGCTGGCT"
    val answer = "AA"
    shortestNonSharedSubstring(text1,text2) should contain (answer)

    val text3 = "GAGCGATGTGTAAAGACGGGCCTAGTGTGTTATTGGTGAGGATCCGCAAATTCCGCTTTAACATCGCTCAACCGTCCACGAGCAGCTCGCGGCTTGTTTGATTTTGCTCACCCGGAACTAAGCCTTCAATAATTGCGGGACACTTCTACTTTCCATGTATTAAAACCTTATAATTCCGAGAGACTGACATTTTAACCCAGGGATTCTCGGGGCTTGTAGGTCACTGGCGTTAGCTCTGCCAATGTTCTTCATCGCCGAAATTAAACTCGTGACCTTGCCCTACGGTTTGAAACGTTGGATTCCTAGCATTCGGTGCAACGGGAGTTCCATACCAGCAGTTAAGGACCGGGTTGCCCCGTCCCACTACGAGCAGCCGTTAGAAAACAGTTCTACCGGAGGCTATCCCGCACCACGGGTTTTCTTAGTGAAAGGGACTGCGCAGCCATCGAAGAGTAGGGGGAGTCAGAGAGAGGCAGGCTTGTTGGGCTGATACATCTAGTTTACTAAATAGCCTTAATGGCGTCCCCCTCTTCGTTGATGCGCGTGGCCTGTGAAATTAGGCAGGGCCCAATGAGCAAGGCTGATTACTATCTAATTGCAGAGCGCAATGCTCTCATATATTATTATCCATGAATCTCATTTCACTAATCAGAAACGTG"
    val text4 = "CACCTCGTCAATACAACAAAAGGCGGCTCGCTTAAAGGGCGCAGCTAGTTCCTCCCCCTCTCATTGGGACATAGTCAACCTGCTAATCCGGATTCGAATGGATTATTCCGTAATTGAACGGTAATTTAGTGAGCTTCGCAGTAAACGATAGATGCGAGCTCTAGCAGGCCACTGACTATATAAACGCCAACACTAGTGCCGTGCATGGACGACTCGATGTACTATAGATTTGCACAGGTATGACCGGAGGAGCGGGACTGCCTAGGCTATAGGGAACGGGGAGTATTGGGAGCCTTTTAGGCCCTCGTCATATCCCTTAACGTTCCCGCGCAGCTAAATTGTGGAACCGGAAAACAATGGATCTGCTTATTTTTGTAGGCTTGGTTAAGCGAAACGGATCAAAATAAACAAAGAATTAATCAATGAACTAACCAACGAAGTAAGCAAGGATATACATAGATTTATTCATTGATCTATCCATCGATGTATGCATGGACACAGACTTACTCACTGACCTACCCACCGACGTACGCACGGAGAGTTAGTCAGTGAGCTAGCCAGCGAGGTAGGCAGGGTTTTCTTTGTTCCTTCGTTGCTTGGTCTCTGTCCCTCCGTCGCTCGGTGTGCCTGCGTGGCTGGGCCCCGCCGGCGCGGGGAAA"
    val answer2 = "AAACC"
    shortestNonSharedSubstring(text3,text4) should contain (answer2)
  }

  it should "recover a string transformed by BWT" in {
    val text = "TTCCTAACG$A"
    val answer = "TACATCACGT$"
    inverseBurrowsWheelerTransform(text) should be (answer)
  }

  it should "match using BWT" in {
    val bwt = "TCCTCTATGAGATCCTATTCTATGAAACCTTCA$GACCAAAATTCTCCGGC"
    val patterns = IndexedSeq[String]("CCT","CAC","GAG","CAG","ATC")
    val answer = IndexedSeq[Int](2,1,1,0,1)
    betterBurrowsWheelerMatching(bwt,patterns) should be (answer)
  }


}
