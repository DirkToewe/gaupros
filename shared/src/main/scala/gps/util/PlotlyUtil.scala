/* This file is part of GauProS.
 *
 * GauProS is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * GauProS is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with GauProS.  If not, see <https://www.gnu.org/licenses/>.
 */

package gps.util

import java.awt.Desktop.getDesktop
import java.nio.file.Files
import java.util.Arrays.asList

/**
  * Created by Dirk Toewe on 27.08.17.
  */
object PlotlyUtil
{
  /** Creates and shows (in the browser) a new <a href="www.plot.ly">Plotly</a> plot.
    *
    * @param layout The plot's layout as json string.
    * @param data The plot's charts as json strings.
    */
  def plot( layout: CharSequence, data: CharSequence* ): Unit =
  {
    val plot = s"""
     |<!DOCTYPE html>
     |<html lang=\"en\">
     |  <head>
     |    <meta charset=\"utf-8\">
     |    <title>Plotly Plot</title>
     |    <script type=\"text/javascript\" src=\"https://cdn.plot.ly/plotly-latest.min.js\"></script>
     |  </head>
     |  <body>
     |    <script type=\"text/javascript\">
     |    'use strict'; {
     |      let div = document.createElement('div');
     |      div.style = 'width: 100%; height: 1024px';
     |      div.innerHTML = 'Creating Plot...';
     |      document.body.appendChild(div);
     |
     |      let
     |        data = [
     |          ${ data mkString ",\n" }
     |        ],
     |        layout = $layout;
     |      div.innerHTML = '';
     |      Plotly.plot(div, data, layout, { showLink: false, modeBarButtonsToRemove: ['sendDataToCloud'] });
     |    }
     |    </script>
     |  </body>
     |</html>
    """.stripMargin

    val tmp = Files.createTempFile("plot_",".html")
    Files.write(tmp, asList(plot) )
    println(tmp)
//    getDesktop.browse(tmp.toUri)
  }
}
