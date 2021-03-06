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
import bioinf.mutations.Mutations
import Mutations._
import bioinf.Answer._
import bioinf.Output._
import scala.io.Source
val it = Source.fromFile("""C:\tmp\dataset_310_2.txt""").getLines()
val text = it.next()
val a = suffixArrayConstruction(text)
val sb = new StringBuilder(a.length * 10)
a.foreach{s =>
  sb.append(s)
  sb.append(", ")
}
sb.delete(sb.length - 2, sb.length)
val answer = sb.mkString
writeStringToFile(answer,"""c:\tmp\results.txt""")
