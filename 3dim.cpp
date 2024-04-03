#include <condition_variable>
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
#include <Magnum/Primitives/Icosphere.h>
#include <Magnum/Shaders/PhongGL.h>
#include <Magnum/Trade/MeshData.h>
#include <mutex>
#include <thread>

#include "arcball/ArcBall.h"
#include "atom.h"

namespace Magnum { namespace Examples {

// Scale changes a particles position in the simulation to the
// particles position in the magnum library.
// VTK has bugs with small scale.  Use scale factor to overcome issues.
// Scale only used to set pseudo position in magnum
const SFloat kScale = 1e-11;    // Number arrived by trial and refinement.
// Real   position = magnum position * scale
// magnum position = real   position / scale


struct SphereInstanceData {
    Matrix4 transformationMatrix;
    Matrix3x3 normalMatrix;
    Color4 color;
};

// Contains logic to display the simulation.
class ThreeDim: public Platform::Application {
public:
  __attribute__((unused)) explicit
  ThreeDim(const Arguments& arguments);

  ~ThreeDim() {  // Destructor
    move_particles_thread_run = false;
    animation_running = false;
    NotifyDrawEvent();  // Inform thread to stop waiting on draw event.
    // std::cout << "\t Requesting MoveParticlesThread() thread to terminate." << std::endl;

    move_particles_thread.join();
  }


protected:
  void viewportEvent(ViewportEvent& event) override;      // Handle window resize
  void keyPressEvent(KeyEvent& event) override;
  void drawEvent() override;                              // Called every frame
  void mousePressEvent(MouseEvent& event) override;
  void mouseReleaseEvent(MouseEvent& event) override;
  void mouseMoveEvent(MouseMoveEvent& event) override;
  void mouseScrollEvent(MouseScrollEvent& event) override;

  void MoveParticlesThread();
  void drawSpheres();

  Containers::Optional<ArcBall> _arcballCamera;
  Matrix4 _projectionMatrix;

  /* Points data as spheres with size */
  Containers::Array<Vector3> _spherePositions;
  Containers::Array<Vector3> _sphereVelocities;
  Float _sphereRadius, _sphereVelocity;
  volatile bool animation_running = true;

  /* Profiling */
  DebugTools::FrameProfilerGL _profiler{
      DebugTools::FrameProfilerGL::Value::FrameTime|
      DebugTools::FrameProfilerGL::Value::CpuDuration, 180};

  /* Spheres rendering */
  GL::Mesh _sphereMesh{NoCreate};
  GL::Buffer _sphereInstanceBuffer{NoCreate};
  Shaders::PhongGL _sphereShader{NoCreate};
  Containers::Array<SphereInstanceData> _sphereInstanceData;

private:
  sew::Atom * atom = new sew::Atom(0);  // Stupid initialization because of compiler error.
  UnsignedInt numSpheres;    // Number of subatomic particles to simulate in the atom.
  std::thread move_particles_thread;
  volatile bool move_particles_thread_run = true;
  void NotifyDrawEvent();
};

using namespace Math::Literals;

ThreeDim::ThreeDim(const Arguments& arguments) : Platform::Application{arguments, NoCreate} {
  Utility::Arguments args;
  args.addOption('s', "spheres", "2")
          .setHelp("spheres", "number of spheres to simulate", "N")
      .addOption('r', "sphere-radius", "0.025")
          .setHelp("sphere-radius", "sphere radius", "R")
      .addOption('v', "sphere-velocity", "0.05")
          .setHelp("sphere-velocity", "sphere velocity", "V")
      .addSkippedPrefix("magnum")
      .parse(arguments.argc, arguments.argv);

  _sphereRadius = args.value<Float>("sphere-radius");
  _sphereVelocity = args.value<Float>("sphere-velocity");
  sew::Atom * atom_ptr;
  // Setup window and parameters
  {
      const Vector2 dpiScaling = this->dpiScaling({});
      Configuration conf;
      conf.setTitle("Magnum Octree Example")
          .setSize({1024+512+256, 1024+128}, dpiScaling)
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
  _spherePositions = Containers::Array<Vector3>{NoInit, numSpheres};
  _sphereVelocities = Containers::Array<Vector3>{NoInit, numSpheres};
  _sphereInstanceData = Containers::Array<SphereInstanceData>{NoInit, numSpheres};
  atom = new sew::Atom(numSpheres);
  atom_ptr = atom->pars[0]->a_;

  _spherePositions[0] = Vector3{
      static_cast<float>(atom->pars[0]->pos[0] / kScale),
      static_cast<float>(atom->pars[0]->pos[1] / kScale),
      static_cast<float>(atom->pars[0]->pos[2] / kScale)};
  for (int i=0; i<numSpheres; i++) {
    sew::Particle *p = atom->pars[i];
    _sphereInstanceData[i].color = Color4{(float)p->color[0]/255,
                                          (float)p->color[1]/255,
                                          (float)p->color[2]/255,
                                          i > numSpheres/2 ? 0.30f : 0.7f};
    for (int j=0; j<3; j++) {
      _spherePositions[i][j] = static_cast<float>(p->pos[j] / kScale);
    }
    float sphere_radius = (p->is_electron || numSpheres <= 4 ? 1 : 1.5) * _sphereRadius;
    _sphereInstanceData[i].transformationMatrix = Matrix4::translation(_spherePositions[i]) *
                                                  Matrix4::scaling(Vector3{sphere_radius}) * 3;
    _sphereInstanceData[i].normalMatrix = _sphereInstanceData[i].transformationMatrix.normalMatrix();
  }
  std::cout << "\t sphere pos: " << _spherePositions[0][0]
              << " electron pos: " << atom->pars[0]->pos[0]
              << "  class atom address " << atom->pars[0]->a_
              << std::endl;
  {        // Rendering spheres / particles.
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
  assert(atom_ptr == atom->pars[0]->a_);
  // Start thread to move particles.
  move_particles_thread = std::thread(&ThreeDim::MoveParticlesThread, this);
  // move_particles_thread.detach(); // Avoids terminal error: terminate called without an active exception
}

std::mutex mutually_exclusive_lock;
std::condition_variable condition_var;

/*
void ThreeDim::MoveParticlesThread() {
  bool just_waited = false;
  while (move_particles_thread_run) {
    while (atom->screen_draw_event_occurred && animation_running) {
      atom->screen_draw_event_occurred = false;
      if (!just_waited)                       // Only used for logging.
        atom->num_drawing_event_already += 1;
      just_waited = false;                    // Only used for logging.

      atom->MoveParticlesThread();
    }
    atom->num_wait_for_drawing_event += 1;     // Only used for logging.
    just_waited = true;                       // Only used for logging.
    {
      // Grab a lock and wait for the next screen draw event.
      std::unique_lock<std::mutex> lk(mutually_exclusive_lock);
      condition_var.wait(lk, [this]{return atom->screen_draw_event_occurred;});
    }
  }
  // std::cout << "\t MoveParticlesThread() thread terminating as requested." << std::endl;
}
*/


// Below runs as a separate thread.  Runs at most once for each frame draw.
void ThreeDim::MoveParticlesThread() {
  while (move_particles_thread_run) {
    // Execute atom move particles on every frame draw event.
    // Don't re-execute atom MoveParticles() without a frame draw event.
    // This enables slowing down the simulation so electrons don't move too fast.
    while (atom->frame_draw_event_occurred && animation_running && move_particles_thread_run) {
      atom->frame_draw_event_occurred = false;
      atom->MoveParticles();
      // Did MoveParticles() executed longer than it took to draw a frame?
      if (atom->frame_draw_event_occurred)
           ++atom->n_times_per_screen_log_MoveParticles_not_compl_before_next_frame_draw_event;
      else ++atom->n_times_per_screen_log_MoveParticles_completed_before_next_frame_draw_event;
    }
    if (!move_particles_thread_run) break;
    {
      // Grab a lock and wait for the next screen draw event.
      std::unique_lock<std::mutex> lk(mutually_exclusive_lock);
      condition_var.wait(lk);
    }
  }
  // std::cout << "\t MoveParticles() thread terminating as requested." << std::endl;
}


void ThreeDim::NotifyDrawEvent() {
  atom->frame_draw_event_occurred = true;
  std::lock_guard<std::mutex> lk(mutually_exclusive_lock);
  condition_var.notify_one();
}

// Draw event is called every frame.
void ThreeDim::drawEvent() {
  NotifyDrawEvent();
  GL::defaultFramebuffer.clear(GL::FramebufferClear::Color|GL::FramebufferClear::Depth);
  _profiler.beginFrame();

  if(animation_running) {
    // atom.MoveParticlesThread();
    for (int i=0; i<numSpheres; i++) {
      for (int j=0; j<3; j++) {
        _spherePositions[i][j] = static_cast<float>(atom->pars[i]->pos[j] / kScale);
      }
    }
  }

  /* Update camera before drawing instances */
  const bool moving = _arcballCamera->updateTransformation();

  drawSpheres();

  _profiler.endFrame();
  // _profiler.printStatistics(300);

  swapBuffers();

  /* If the camera is moving or the animation is running, redraw immediately */
  if(moving || animation_running) redraw();
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
    if(event.key() == KeyEvent::Key::D) {
        atom->DtLoggingToggle();
    } else if(event.key() == KeyEvent::Key::E) {
        atom->EnergyLoggingToggle();
    } else if(event.key() == KeyEvent::Key::F) {
        atom->FastLoggingToggle();
    } else if(event.key() == KeyEvent::Key::I) {
        atom->IterationsLoggingToggle();
    } else if(event.key() == KeyEvent::Key::O) {
        if(_profiler.isEnabled()) _profiler.enable();
    } else if(event.key() == KeyEvent::Key::P) {
        atom->PositionLoggingToggle();
    } else if(event.key() == KeyEvent::Key::Q) {
        if(_profiler.isEnabled()) _profiler.disable();
        else _profiler.enable();
    } else if(event.key() == KeyEvent::Key::R) {
        _arcballCamera->reset();
    } else if(event.key() == KeyEvent::Key::V) {
        atom->VelocityLoggingToggle();
    } else if(event.key() == KeyEvent::Key::Y) {
      atom->PercentEnergyDissipatedLoggingToggle();
    } else if(event.key() == KeyEvent::Key::Z) {
        atom->FrameDrawStatisticsLoggingToggle();
    } else if(event.key() == KeyEvent::Key::Space) {
        animation_running ^= true;

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
