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

import org.scalajs.jsenv.nodejs.NodeJSEnv
import sbtcrossproject.CrossPlugin.autoImport.{crossProject, CrossType}

name := "gauprosRoot"

// TO TEST RUN
//   - gauprosJVM/test
//   - gauprosJS/testint

//lazy val gauprosBuild = project.in( file(".") )
//  .aggregate(gauprosJS, gauprosJVM)
//  .settings()


lazy val gaupros = crossProject(JSPlatform, JVMPlatform)
  .crossType(CrossType.Full)
  .in( file(".") )
  .settings(
     name := "gaupros",
     version := "0.3.0",
     scalaVersion := "2.12.8",
     scalacOptions ++= Seq(
      "-feature",
      "-deprecation"
     ),

     testFrameworks += new TestFramework("utest.runner.Framework"),

     libraryDependencies ++= Seq(
       "com.lihaoyi" %%% "utest" % "0.6.6" % "test"
     )
   )
  .jsSettings(
     jsEnv := new NodeJSEnv( NodeJSEnv.Config().withArgs("--max_old_space_size=3584" :: Nil) ),
     scalaJSUseMainModuleInitializer := true,
     scalaJSStage := FullOptStage
   )
  .jvmSettings(
   )


lazy val gauprosJS = gaupros.js
lazy val gauprosJVM= gaupros.jvm
