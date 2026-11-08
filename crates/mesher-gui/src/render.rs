//! The 3D viewport: a wgpu render pass into an offscreen colour and depth
//! target, shown in egui as an image. Flat-shaded triangles with a headlight,
//! lines for the wireframe and the defect overlay, and an optional clipping
//! plane.

use eframe::egui;
use eframe::egui_wgpu::{self, wgpu};
use wgpu::util::DeviceExt;

#[repr(C)]
#[derive(Debug, Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Vertex {
    pub pos: [f32; 3],
    pub normal: [f32; 3],
    /// Linear RGBA.
    pub color: [f32; 4],
}

#[repr(C)]
#[derive(Debug, Clone, Copy, bytemuck::Pod, bytemuck::Zeroable)]
pub struct Uniforms {
    pub view_proj: [[f32; 4]; 4],
    /// Direction towards the light (xyz).
    pub light: [f32; 4],
    /// Clipping plane: points with `dot(xyz, p) > w` are cut away.
    pub clip: [f32; 4],
    /// x: clipping on (1) or off (0).
    pub flags: [f32; 4],
}

const COLOR: wgpu::TextureFormat = wgpu::TextureFormat::Rgba8UnormSrgb;
const DEPTH: wgpu::TextureFormat = wgpu::TextureFormat::Depth32Float;

const SHADER: &str = r"
struct U { view_proj: mat4x4<f32>, light: vec4<f32>, clip: vec4<f32>, flags: vec4<f32> };
@group(0) @binding(0) var<uniform> u: U;

struct VIn {
    @location(0) pos: vec3<f32>,
    @location(1) normal: vec3<f32>,
    @location(2) color: vec4<f32>,
};
struct VOut {
    @builtin(position) clip_pos: vec4<f32>,
    @location(0) normal: vec3<f32>,
    @location(1) color: vec4<f32>,
    @location(2) world: vec3<f32>,
};

@vertex
fn vs(v: VIn) -> VOut {
    var o: VOut;
    o.clip_pos = u.view_proj * vec4<f32>(v.pos, 1.0);
    o.normal = v.normal;
    o.color = v.color;
    o.world = v.pos;
    return o;
}

fn clipped(world: vec3<f32>) -> bool {
    return u.flags.x > 0.5 && dot(u.clip.xyz, world) > u.clip.w;
}

@fragment
fn fs_tri(i: VOut, @builtin(front_facing) front: bool) -> @location(0) vec4<f32> {
    if (clipped(i.world)) { discard; }
    var n = normalize(i.normal);
    if (!front) { n = -n; }
    let l = abs(dot(n, normalize(u.light.xyz)));
    // Back faces (seen through a clip) are darker.
    var shade = 0.3 + 0.7 * l;
    if (!front) { shade = shade * 0.55; }
    return vec4<f32>(i.color.rgb * shade, 1.0);
}

@fragment
fn fs_line(i: VOut) -> @location(0) vec4<f32> {
    if (clipped(i.world)) { discard; }
    return i.color;
}
";

pub struct Viewport {
    tris: wgpu::RenderPipeline,
    lines: wgpu::RenderPipeline,
    uniforms: wgpu::Buffer,
    bind_group: wgpu::BindGroup,
    tri_buffer: Option<(wgpu::Buffer, u32)>,
    line_buffer: Option<(wgpu::Buffer, u32)>,
    target: Option<Target>,
}

struct Target {
    size: [u32; 2],
    color: wgpu::TextureView,
    depth: wgpu::TextureView,
    id: egui::TextureId,
}

impl Viewport {
    #[must_use]
    pub fn new(rs: &egui_wgpu::RenderState) -> Viewport {
        let device = &rs.device;
        let shader = device.create_shader_module(wgpu::ShaderModuleDescriptor {
            label: Some("viewport"),
            source: wgpu::ShaderSource::Wgsl(SHADER.into()),
        });
        let uniforms = device.create_buffer(&wgpu::BufferDescriptor {
            label: Some("viewport uniforms"),
            size: std::mem::size_of::<Uniforms>() as u64,
            usage: wgpu::BufferUsages::UNIFORM | wgpu::BufferUsages::COPY_DST,
            mapped_at_creation: false,
        });
        let layout = device.create_bind_group_layout(&wgpu::BindGroupLayoutDescriptor {
            label: Some("viewport"),
            entries: &[wgpu::BindGroupLayoutEntry {
                binding: 0,
                visibility: wgpu::ShaderStages::VERTEX_FRAGMENT,
                ty: wgpu::BindingType::Buffer {
                    ty: wgpu::BufferBindingType::Uniform,
                    has_dynamic_offset: false,
                    min_binding_size: None,
                },
                count: None,
            }],
        });
        let bind_group = device.create_bind_group(&wgpu::BindGroupDescriptor {
            label: Some("viewport"),
            layout: &layout,
            entries: &[wgpu::BindGroupEntry {
                binding: 0,
                resource: uniforms.as_entire_binding(),
            }],
        });
        let pipeline_layout = device.create_pipeline_layout(&wgpu::PipelineLayoutDescriptor {
            label: Some("viewport"),
            bind_group_layouts: &[Some(&layout)],
            immediate_size: 0,
        });
        let attributes = wgpu::vertex_attr_array![0 => Float32x3, 1 => Float32x3, 2 => Float32x4];
        let buffers = [Some(wgpu::VertexBufferLayout {
            array_stride: std::mem::size_of::<Vertex>() as u64,
            step_mode: wgpu::VertexStepMode::Vertex,
            attributes: &attributes,
        })];
        let pipeline = |topology, fragment: &str, bias: i32| {
            device.create_render_pipeline(&wgpu::RenderPipelineDescriptor {
                label: Some("viewport"),
                layout: Some(&pipeline_layout),
                vertex: wgpu::VertexState {
                    module: &shader,
                    entry_point: Some("vs"),
                    compilation_options: wgpu::PipelineCompilationOptions::default(),
                    buffers: &buffers,
                },
                primitive: wgpu::PrimitiveState {
                    topology,
                    ..wgpu::PrimitiveState::default()
                },
                depth_stencil: Some(wgpu::DepthStencilState {
                    format: DEPTH,
                    depth_write_enabled: Some(true),
                    depth_compare: Some(wgpu::CompareFunction::LessEqual),
                    stencil: wgpu::StencilState::default(),
                    // Push faces back a little so their edges draw on top.
                    bias: wgpu::DepthBiasState {
                        constant: bias,
                        slope_scale: if bias > 0 { 1.0 } else { 0.0 },
                        clamp: 0.0,
                    },
                }),
                multisample: wgpu::MultisampleState::default(),
                fragment: Some(wgpu::FragmentState {
                    module: &shader,
                    entry_point: Some(fragment),
                    compilation_options: wgpu::PipelineCompilationOptions::default(),
                    targets: &[Some(wgpu::ColorTargetState {
                        format: COLOR,
                        blend: Some(wgpu::BlendState::ALPHA_BLENDING),
                        write_mask: wgpu::ColorWrites::ALL,
                    })],
                }),
                multiview_mask: None,
                cache: None,
            })
        };
        Viewport {
            tris: pipeline(wgpu::PrimitiveTopology::TriangleList, "fs_tri", 2),
            lines: pipeline(wgpu::PrimitiveTopology::LineList, "fs_line", 0),
            uniforms,
            bind_group,
            tri_buffer: None,
            line_buffer: None,
            target: None,
        }
    }

    /// Replace the geometry: triangle-list and line-list vertices.
    pub fn set_geometry(&mut self, rs: &egui_wgpu::RenderState, tris: &[Vertex], lines: &[Vertex]) {
        let make = |data: &[Vertex]| {
            (!data.is_empty()).then(|| {
                let buffer = rs
                    .device
                    .create_buffer_init(&wgpu::util::BufferInitDescriptor {
                        label: Some("viewport geometry"),
                        contents: bytemuck::cast_slice(data),
                        usage: wgpu::BufferUsages::VERTEX,
                    });
                (buffer, data.len() as u32)
            })
        };
        self.tri_buffer = make(tris);
        self.line_buffer = make(lines);
    }

    /// Render a frame of `size` pixels and return the egui texture holding it.
    pub fn render(
        &mut self,
        rs: &egui_wgpu::RenderState,
        size: [u32; 2],
        uniforms: &Uniforms,
        background: [f64; 4],
    ) -> egui::TextureId {
        let size = [size[0].max(1), size[1].max(1)];
        let device = &rs.device;
        if self.target.as_ref().is_none_or(|t| t.size != size) {
            let texture = |format, usage| {
                device
                    .create_texture(&wgpu::TextureDescriptor {
                        label: Some("viewport target"),
                        size: wgpu::Extent3d {
                            width: size[0],
                            height: size[1],
                            depth_or_array_layers: 1,
                        },
                        mip_level_count: 1,
                        sample_count: 1,
                        dimension: wgpu::TextureDimension::D2,
                        format,
                        usage,
                        view_formats: &[],
                    })
                    .create_view(&wgpu::TextureViewDescriptor::default())
            };
            let color = texture(
                COLOR,
                wgpu::TextureUsages::RENDER_ATTACHMENT | wgpu::TextureUsages::TEXTURE_BINDING,
            );
            let depth = texture(DEPTH, wgpu::TextureUsages::RENDER_ATTACHMENT);
            let mut renderer = rs.renderer.write();
            let id = match self.target.take() {
                Some(old) => {
                    renderer.update_egui_texture_from_wgpu_texture(
                        device,
                        &color,
                        wgpu::FilterMode::Linear,
                        old.id,
                    );
                    old.id
                }
                None => renderer.register_native_texture(device, &color, wgpu::FilterMode::Linear),
            };
            self.target = Some(Target {
                size,
                color,
                depth,
                id,
            });
        }
        let Some(target) = &self.target else {
            return egui::TextureId::default();
        };
        rs.queue
            .write_buffer(&self.uniforms, 0, bytemuck::bytes_of(uniforms));
        let mut encoder = device.create_command_encoder(&wgpu::CommandEncoderDescriptor {
            label: Some("viewport"),
        });
        {
            let [r, g, b, a] = background;
            let mut pass = encoder.begin_render_pass(&wgpu::RenderPassDescriptor {
                label: Some("viewport"),
                color_attachments: &[Some(wgpu::RenderPassColorAttachment {
                    view: &target.color,
                    depth_slice: None,
                    resolve_target: None,
                    ops: wgpu::Operations {
                        load: wgpu::LoadOp::Clear(wgpu::Color { r, g, b, a }),
                        store: wgpu::StoreOp::Store,
                    },
                })],
                depth_stencil_attachment: Some(wgpu::RenderPassDepthStencilAttachment {
                    view: &target.depth,
                    depth_ops: Some(wgpu::Operations {
                        load: wgpu::LoadOp::Clear(1.0),
                        store: wgpu::StoreOp::Store,
                    }),
                    stencil_ops: None,
                }),
                timestamp_writes: None,
                occlusion_query_set: None,
                multiview_mask: None,
            });
            pass.set_bind_group(0, &self.bind_group, &[]);
            if let Some((buffer, count)) = &self.tri_buffer {
                pass.set_pipeline(&self.tris);
                pass.set_vertex_buffer(0, buffer.slice(..));
                pass.draw(0..*count, 0..1);
            }
            if let Some((buffer, count)) = &self.line_buffer {
                pass.set_pipeline(&self.lines);
                pass.set_vertex_buffer(0, buffer.slice(..));
                pass.draw(0..*count, 0..1);
            }
        }
        rs.queue.submit([encoder.finish()]);
        target.id
    }
}
