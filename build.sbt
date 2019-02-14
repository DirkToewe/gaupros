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
     version := "0.2.0",
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
