import org.scalajs.jsenv.nodejs.NodeJSEnv
// import sbtcrossproject.{CrossType, crossProject}

// enablePlugins(ScalaJSPlugin)

name := "GauProS"
scalaVersion in ThisBuild := "2.12.4"

scalaJSStage in Global := FullOptStage
// scalaJSStage in Global := FastOptStage
scalacOptions in Global ++= Seq("-feature", "-deprecation")

lazy val root = project.in( file(".") )
  .aggregate(gauprosJS, gauprosJVM)
  .settings()

lazy val gaupros = crossProject //JSPlatform, JVMPlatform)
  .in( file(".") )
//  .crossType(CrossType.Pure)
  .settings(
      name := "GauProS",
      version := "0.2",
      scalaVersion := "2.12.4",

      testFrameworks += new TestFramework("utest.runner.Framework"),

      libraryDependencies += "com.lihaoyi" %%% "utest" % "0.6.2" % "test"
  )
  .jsSettings(
      jsEnv := new NodeJSEnv( NodeJSEnv.Config().withArgs("--max_old_space_size=12288" :: Nil) ),
      scalaJSUseMainModuleInitializer := true
  )
  .jvmSettings(
  )
  .in( file(".") )

lazy val gauprosJS = gaupros.js
lazy val gauprosJVM= gaupros.jvm
