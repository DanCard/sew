#include <Corrade/Containers/Optional.h>
#include <Corrade/Containers/GrowableArray.h>
#include <Corrade/Containers/Pointer.h>
#include <Corrade/Utility/Arguments.h>
#include <Magnum/DebugTools/FrameProfiler.h>
#include <Magnum/GL/DefaultFramebuffer.h>
#include <Magnum/GL/Renderer.h>
#include <Magnum/GL/Mesh.h>
#include <Magnum/Math/Color.h>
#include <Magnum/MeshTools/Compile.h>
#include <Magnum/Platform/Sdl2Application.h>
#include <Magnum/Primitives/Cube.h>
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Trade/MeshData.h>

#include "arcball/ArcBall.h"
#include "particles.cpp"

namespace Magnum { namespace Examples {

const bool useStandingEMWave = true;
// Scale changes a particles position in the simulation to the
// particles position in the magnum library.
// VTK has bugs with small scale.  Use scale factor to overcome issues.
// Scale only used to set pseudo position in magnum
const double kScale = 1e-11;    // Number arrived by trial and refinement.
// Real   position = magnum position * scale
// magnum position = real   position / scale


struct SphereInstanceData {
    Matrix4 transformationMatrix;
    Matrix3x3 normalMatrix;
    Color4 color;
};

class ThreeDim: public Platform::Application {
    public:
    __attribute__((unused)) explicit
    ThreeDim(const Arguments& arguments);

    protected:
        Particles atom = Particles(0);  // Stupid initialization because of compiler error.
        UnsignedInt numSpheres;  // Number of subatomic particles to simulate in the atom.
        void viewportEvent(ViewportEvent& event) override;      // Handle window resize
        void keyPressEvent(KeyEvent& event) override;
        void drawEvent() override;                              // Called every frame
        void mousePressEvent(MouseEvent& event) override;
        void mouseReleaseEvent(MouseEvent& event) override;
        void mouseMoveEvent(MouseMoveEvent& event) override;
        void mouseScrollEvent(MouseScrollEvent& event) override;

        void moveParticles();
        void drawSpheres();

        Containers::Optional<ArcBall> _arcballCamera;
        Matrix4 _projectionMatrix;

        /* Points data as spheres with size */
        Containers::Array<Vector3> _spherePositions;
        Containers::Array<Vector3> _sphereVelocities;
        Float _sphereRadius, _sphereVelocity;
        bool _animation = true;

        /* Profiling */
        DebugTools::FrameProfilerGL _profiler{
            DebugTools::FrameProfilerGL::Value::FrameTime|
            DebugTools::FrameProfilerGL::Value::CpuDuration, 180};

        /* Spheres rendering */
        GL::Mesh _sphereMesh{NoCreate};
        GL::Buffer _sphereInstanceBuffer{NoCreate};
        Shaders::PhongGL _sphereShader{NoCreate};
        Containers::Array<SphereInstanceData> _sphereInstanceData;
};

using namespace Math::Literals;

ThreeDim::ThreeDim(const Arguments& arguments) : Platform::Application{arguments, NoCreate} {
  Utility::Arguments args;
  args.addOption('s', "spheres", "4")
          .setHelp("spheres", "number of spheres to simulate", "N")
      .addOption('r', "sphere-radius", "0.025")
          .setHelp("sphere-radius", "sphere radius", "R")
      .addOption('v', "sphere-velocity", "0.05")
          .setHelp("sphere-velocity", "sphere velocity", "V")
      .addSkippedPrefix("magnum")
      .parse(arguments.argc, arguments.argv);

  _sphereRadius = args.value<Float>("sphere-radius");
  _sphereVelocity = args.value<Float>("sphere-velocity");

  // Setup window and parameters
  {
      const Vector2 dpiScaling = this->dpiScaling({});
      Configuration conf;
      conf.setTitle("Magnum Octree Example")
          .setSize({1024+512, 1024}, dpiScaling)
          .setWindowFlags(Configuration::WindowFlag::Resizable);
      Debug{} << "size:" << conf.size() << "dpiScaling:" << dpiScaling << "max:" << dpiScaling.max();
      GLConfiguration glConf;
      glConf.setSampleCount(dpiScaling.max() < 2.0f ? 8 : 2);
      // glConf.set
      if(!tryCreate(conf, glConf)) {
          create(conf, glConf.setSampleCount(0));
      }

      SDL_SetWindowPosition(this->window(), 40, 0);

      GL::Renderer::enable(GL::Renderer::Feature::DepthTest);
      // https://doc.magnum.graphics/magnum/classMagnum_1_1Shaders_1_1PhongGL.html#Shaders-PhongGL-alpha
      GL::Renderer::enable(GL::Renderer::Feature::Blending);
      GL::Renderer::enable(GL::Renderer::Feature::FaceCulling);
      GL::Renderer::setBlendFunction(
          GL::Renderer::BlendFunction::SourceAlpha,
          GL::Renderer::BlendFunction::OneMinusSourceAlpha);

      /* Loop at 60 Hz max */
      setSwapInterval(1);
      setMinimalLoopPeriod(16);   // 16 milliseconds.  60 Hz = 16.6667 milliseconds
  }

  // Setup camera
  {
      const Vector3 eye = Vector3::zAxis(5.0f);
      const Vector3 viewCenter;
      const Vector3 up = Vector3::yAxis();
      const Deg fov = 45.0_degf;
      _arcballCamera.emplace(eye, viewCenter, up, fov, windowSize());
      _arcballCamera->setLagging(0.85f);

      _projectionMatrix = Matrix4::perspectiveProjection(fov,
          Vector2{framebufferSize()}.aspectRatio(), 0.01f, 100.0f);
  }

  numSpheres = args.value<UnsignedInt>("spheres");
  const Vector3 tmpPos = Vector3(std::rand(), std::rand(), std::rand()) / Float(RAND_MAX);
  _spherePositions = Containers::Array<Vector3>{NoInit, numSpheres};
  _sphereVelocities = Containers::Array<Vector3>{NoInit, numSpheres};
  _sphereInstanceData = Containers::Array<SphereInstanceData>{NoInit, numSpheres};
  // if (useStandingEMWave) {
      atom = Particles(numSpheres);
      _spherePositions[0] = Vector3{
          static_cast<float>(atom.pars[0]->pos[0] / kScale),
          static_cast<float>(atom.pars[0]->pos[1] / kScale),
          static_cast<float>(atom.pars[0]->pos[2] / kScale)};
      std::cout << "\t sphere pos: " << _spherePositions[0][0]
                  << " electron pos: " << atom.pars[0]->pos[0]
                  << std::endl;
      _spherePositions[1] = Vector3{0.0f};
      _sphereInstanceData[0].color = Color4{1.0F, 0.0f, 0.0f, 0.5f};
      _sphereInstanceData[1].color = Color4{0.0f, 0.0f, 1.0f, 0.5f};
      /*
  } else {
      // Setup points (render as spheres)

      for(std::size_t i = 0; i < numSpheres; ++i) {
          const Vector3 tmpVel = Vector3(std::rand(), std::rand(), std::rand())/
              Float(RAND_MAX);
          _spherePositions[i] = tmpPos*2.0f - Vector3{1.0f};
          _spherePositions[i].y() *= 0.5f;
          _sphereVelocities[i] = (tmpVel*2.0f - Vector3{1.0f}).resized(_sphereVelocity);
          _sphereInstanceData[i].color = Color3{tmpPos};  // Should that be Color4?
      }
  }
      */
  {        // Rendering spheres / particles.
      for(std::size_t i = 0; i < numSpheres; ++i) {
          /* Fill in the instance data. Most of these stays the same, except
              for the translation */
          // If we don't multiply by the _sphereRadius, the spheres will be too big.
          // No idea why.
          _sphereInstanceData[i].transformationMatrix =
              Matrix4::translation(_spherePositions[i]) *
              Matrix4::scaling(Vector3{_sphereRadius}) * 3;
          _sphereInstanceData[i].normalMatrix =
              _sphereInstanceData[i].transformationMatrix.normalMatrix();
      }

      _sphereShader = Shaders::PhongGL{Shaders::PhongGL::Configuration{}
          .setFlags(Shaders::PhongGL::Flag::VertexColor|
                      Shaders::PhongGL::Flag::InstancedTransformation)};
      _sphereInstanceBuffer = GL::Buffer{};
      _sphereMesh = MeshTools::compile(Primitives::icosphereSolid(2));
      _sphereMesh.addVertexBufferInstanced(_sphereInstanceBuffer, 1, 0,
          Shaders::PhongGL::TransformationMatrix{},
          Shaders::PhongGL::NormalMatrix{},
  //      Shaders::PhongGL::Color3{});      // Should that be Color4?
          Shaders::PhongGL::Color4{});      // Should that be Color4?
      _sphereMesh.setInstanceCount(_sphereInstanceData.size());
  }
}

    // Draw event is called every frame.
void ThreeDim::drawEvent() {
  GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);
  _profiler.beginFrame();

  if(_animation) {
    if (useStandingEMWave) {
      atom.moveParticles();
      for (int i=0; i<numSpheres; i++) {
        for (int j=0; j<3; j++) {
          _spherePositions[i][j] = static_cast<float>(atom.pars[i]->pos[j] / kScale);
        }
      }
    } else {
        moveParticles();
    }
  }

  /* Update camera before drawing instances */
  const bool moving = _arcballCamera->updateTransformation();

  drawSpheres();

  _profiler.endFrame();
  // _profiler.printStatistics(300);

  swapBuffers();

  /* If the camera is moving or the animation is running, redraw immediately */
  if(moving || _animation) redraw();
}

void ThreeDim::moveParticles() {
    constexpr Float dt = 1.0f/120.0f;

    for(std::size_t i = 0; i < _spherePositions.size(); ++i) {
        Vector3 pos = _spherePositions[i] + _sphereVelocities[i] * dt;
        for(std::size_t j = 0; j < 3; ++j) {
            if(pos[j] < -1.0f || pos[j] > 1.0f)
                _sphereVelocities[i][j] = -_sphereVelocities[i][j];
            pos[j] = Math::clamp(pos[j], -1.0f, 1.0f);
        }

        _spherePositions[i] = pos;
    }
}

void ThreeDim::drawSpheres() {
    // Loop through all the spheres and update their transformation matrix
    for(std::size_t i = 0; i != _spherePositions.size(); ++i)
        _sphereInstanceData[i].transformationMatrix.translation() =
            _spherePositions[i];

    _sphereInstanceBuffer.setData(_sphereInstanceData, GL::BufferUsage::DynamicDraw);
    _sphereShader
        .setProjectionMatrix(_projectionMatrix)
        .setTransformationMatrix(_arcballCamera->viewMatrix())
        .setNormalMatrix(_arcballCamera->viewMatrix().normalMatrix())
        .draw(_sphereMesh);
}

// Handle window resize.
void ThreeDim::viewportEvent(ViewportEvent& event) {
    GL::defaultFramebuffer.setViewport({{}, event.framebufferSize()});
    _arcballCamera->reshape(event.windowSize());

    _projectionMatrix = Matrix4::perspectiveProjection(_arcballCamera->fov(),
        Vector2{event.framebufferSize()}.aspectRatio(), 0.01f, 100.0f);
}

void ThreeDim::keyPressEvent(KeyEvent& event) {
    if(event.key() == KeyEvent::Key::O) {
        if(_profiler.isEnabled()) _profiler.enable();

    } else if(event.key() == KeyEvent::Key::P) {
        if(_profiler.isEnabled()) _profiler.disable();
        else _profiler.enable();

    } else if(event.key() == KeyEvent::Key::R) {
        _arcballCamera->reset();

    } else if(event.key() == KeyEvent::Key::Space) {
        _animation ^= true;

    } else return;

    event.setAccepted();
    redraw();
}

void ThreeDim::mousePressEvent(MouseEvent& event) {
    /* Enable mouse capture so the mouse can drag outside of the window */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_TRUE);
    _arcballCamera->initTransformation(event.position());
    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void ThreeDim::mouseReleaseEvent(MouseEvent&) {
    /* Disable mouse capture again */
    /** @todo replace once https://github.com/mosra/magnum/pull/419 is in */
    SDL_CaptureMouse(SDL_FALSE);
}

void ThreeDim::mouseMoveEvent(MouseMoveEvent& event) {
    if(!event.buttons()) return;

    if(event.modifiers() & MouseMoveEvent::Modifier::Shift)
        _arcballCamera->translate(event.position());
    else _arcballCamera->rotate(event.position());

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

void ThreeDim::mouseScrollEvent(MouseScrollEvent& event) {
    const Float delta = event.offset().y();
    if(Math::abs(delta) < 1.0e-2f) return;

    _arcballCamera->zoom(delta);

    event.setAccepted();
    redraw(); /* camera has changed, redraw! */
}

}}

MAGNUM_APPLICATION_MAIN(Magnum::Examples::ThreeDim)